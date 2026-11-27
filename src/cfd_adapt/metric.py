"""Double P1 recovery, spectral regularisation and discrete metric construction."""

from dataclasses import dataclass
from numbers import Integral

import numpy as np

from .mesh import Mesh


@dataclass(frozen=True)
class Recovery:
    element_gradient: np.ndarray
    nodal_gradient: np.ndarray
    raw_hessian: np.ndarray
    hessian: np.ndarray


@dataclass(frozen=True)
class MetricParameters:
    epsilon: float
    h_min: float
    h_max: float
    max_aspect: float
    p: int = 1
    target_nodes: int | None = None
    target_complexity: float | None = None
    complexity_rtol: float = 1e-10
    bisection_steps: int = 100

    def __post_init__(self):
        scalars = (self.epsilon, self.h_min, self.h_max, self.max_aspect, self.complexity_rtol)
        if not np.isfinite(scalars).all() or not (
            self.epsilon > 0 and 0 < self.h_min <= self.h_max and self.max_aspect >= 1
            and 0 < self.complexity_rtol < 1
        ):
            raise ValueError("invalid curvature floor, metric limits or complexity tolerance")
        for label, value in (("p", self.p), ("bisection_steps", self.bisection_steps)):
            if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
                raise ValueError(f"{label} must be a positive integer")
        if (self.target_nodes is None) == (self.target_complexity is None):
            raise ValueError("select exactly one target: nodes or complexity")
        if self.target_nodes is not None and (
            isinstance(self.target_nodes, bool) or not isinstance(self.target_nodes, Integral)
            or self.target_nodes < 1
        ):
            raise ValueError("target_nodes must be a positive integer")
        if not np.isfinite(self.complexity) or self.complexity <= 0:
            raise ValueError("target complexity must be positive and finite")

    @property
    def complexity(self):
        """MA-C01: approximate conversion from the requested node count."""
        if self.target_nodes is None:
            return float(self.target_complexity)
        return (np.sqrt(3) / 2) * self.target_nodes


@dataclass(frozen=True)
class MetricResult:
    hessian_eigenvalues: np.ndarray
    principal_directions: np.ndarray
    curvature_eigenvalues: np.ndarray
    curvature: np.ndarray
    optimal: np.ndarray
    bounded: np.ndarray
    optimal_eigenvalues: np.ndarray
    bounded_eigenvalues: np.ndarray
    optimal_complexity: float
    bounded_complexity: float
    feasible_complexities: tuple[float, float]
    log_global_scaling: float
    log_beta: float
    bisection_iterations: int

    @property
    def target_sizes(self):
        return 1 / np.sqrt(self.bounded_eigenvalues)

    @property
    def metric_density(self):
        return np.exp(0.5 * np.log(self.bounded_eigenvalues).sum(axis=1))

    @property
    def anisotropy(self):
        return self.target_sizes[:, 0] / self.target_sizes[:, 1]


def element_gradient(mesh: Mesh, values):
    """MA-R01: solve the P1 edge-difference system.

    Input (N,) or (N, fields). Output (T, 2) or (T, fields, 2).
    Each nondegenerate triangle supplies two physical coordinate differences.
    """
    values = np.asarray(values, dtype=float)
    if values.ndim not in (1, 2) or len(values) != len(mesh.coordinates) or not np.isfinite(values).all():
        raise ValueError("expected finite nodal values with shape (N,) or (N, fields)")
    scalar = values.ndim == 1
    if scalar:
        values = values[:, None]
    points = mesh.coordinates[mesh.triangles]
    local_values = values[mesh.triangles]
    gradient = np.linalg.solve(points[:, 1:] - points[:, :1], local_values[:, 1:] - local_values[:, :1])
    gradient = gradient.transpose(0, 2, 1)
    return gradient[:, 0] if scalar else gradient


def nodal_average(mesh: Mesh, values):
    """MA-R02 and eq:lumped_nodal_area: area-weight (T, ...) to (N, ...)."""
    values = np.asarray(values, dtype=float)
    if values.ndim == 0 or len(values) != len(mesh.triangles) or not np.isfinite(values).all():
        raise ValueError("expected finite element values with shape (T, ...)")
    dimensions = (1,) * (values.ndim - 1)
    weighted = values * mesh.triangle_areas.reshape(-1, *dimensions)
    result = np.zeros((len(mesh.coordinates), *values.shape[1:]))
    for nodes in mesh.triangles.T:
        np.add.at(result, nodes, weighted)
    return result / (3 * mesh.nodal_areas.reshape(-1, *dimensions))


def recover_hessian(mesh: Mesh, density):
    """MA-R03: P1 differentiation and area recovery, followed by symmetrisation.

    Density is (N,), with element/nodal gradients (T,2)/(N,2)
    and nodal Jacobian and Hessian (N,2,2), in density per length squared.
    """
    if np.shape(density) != (len(mesh.coordinates),):
        raise ValueError("density must have one scalar per vertex")
    g_t = element_gradient(mesh, density)
    g_i = nodal_average(mesh, g_t)
    d_i = nodal_average(mesh, element_gradient(mesh, g_i))
    return Recovery(g_t, g_i, d_i, 0.5 * d_i + 0.5 * d_i.swapaxes(1, 2))


def regularise_curvature(hessian, epsilon):
    """eq: regCurv: signed eigenvalues, basis and floored magnitudes per vertex.

    Spectrum (N,2) and basis (N,2,2) are ordered by increasing floored curvature.
    """
    hessian = np.asarray(hessian, dtype=float)
    if (hessian.ndim != 3 or hessian.shape[1:] != (2, 2) or not np.isfinite(hessian).all()
            or not np.array_equal(hessian, hessian.swapaxes(1, 2))):
        raise ValueError("expected finite symmetric Hessians of shape (N, 2, 2)")
    if not np.isfinite(epsilon) or epsilon <= 0:
        raise ValueError("epsilon must be positive and finite")
    eigenvalues, basis = np.linalg.eigh(hessian)
    curvature = np.maximum(np.abs(eigenvalues), epsilon)
    order = np.argsort(curvature, axis=1, kind="stable")
    return (
        np.take_along_axis(eigenvalues, order, axis=1),
        np.take_along_axis(basis, order[:, None, :], axis=2),
        np.take_along_axis(curvature, order, axis=1),
    )


def tensor(basis, spectrum):
    result = (basis * spectrum[:, None, :]) @ basis.swapaxes(1, 2)
    return 0.5 * result + 0.5 * result.swapaxes(1, 2)


def logsumexp(values):
    maximum = np.max(values)
    return maximum + np.log(np.exp(values - maximum).sum())


def optimal_log_spectrum(curvature, areas, p, target):
    """eq:discrete_global_scaling, eq:discrete_optimal_metric.

    Inputs are positive curvature eigenvalues (N,2), nodal areas (N,),
    positive integer p and dimensionless target. Returns log spectrum and log s.
    Using log(det K) = log(k_1) + log(k_2) avoids determinant overflow.
    """
    log_k = np.log(curvature)
    log_det_k = log_k.sum(axis=1)
    half_power = 0.5 / (p + 1)
    log_s = np.log(target) - logsumexp(np.log(areas) + (0.5 - half_power) * log_det_k)
    return log_s - half_power * log_det_k[:, None] + log_k, float(log_s)


def bounded_log_spectrum(log_spectrum, log_beta, log_min, log_max, log_aspect):
    """eq:bounded_metric_beta, eq:clipped_metric_eigenvalues, eq:bounded_aspect_ratio.

    Scale (N,2) ascending spectrum, clip both values, then lift the smaller.
    Logarithms express multiplication and the aspect-ratio bound as additions.
    """
    bounded = np.clip(log_spectrum + log_beta, log_min, log_max)
    bounded[:, 0] = np.maximum(bounded[:, 0], bounded[:, 1] - 2 * log_aspect)
    return bounded


def recover_complexity(log_spectrum, areas, parameters: MetricParameters):
    """eq:bounded_metric_complexity, eq:feasible_complexity_interval.

    Bisect log(beta) between full saturation bounds. Return the bounded log
    spectrum (N,2), log(beta), complexity, feasibility interval and iteration count.
    """
    log_min = -2 * np.log(parameters.h_max)
    log_max = -2 * np.log(parameters.h_min)
    log_aspect = np.log(parameters.max_aspect)
    log_areas = np.log(areas)
    target = parameters.complexity
    log_target = np.log(target)
    log_total_area = np.log(areas.sum())
    limits = (log_total_area + log_min, log_total_area + log_max)
    rtol = parameters.complexity_rtol
    lower_error, upper_error = (bound - log_target for bound in limits)
    # The feasibility allowance addresses floating-point evaluation of the endpoints.
    slack = 8 * np.finfo(float).eps * max(1, abs(log_target), *map(abs, limits))
    with np.errstate(over="ignore", under="ignore"):
        interval = tuple(float(np.exp(bound)) for bound in limits)
        represented_bounds = np.exp([log_min, log_max])
    if np.any(represented_bounds <= 0) or not np.isfinite(represented_bounds).all():
        raise ValueError("metric size bounds exceed the representable floating-point spectrum")
    if lower_error > slack or upper_error < -slack:
        raise ValueError(f"infeasible target complexity {target:g}; feasible interval {interval}")
    low = np.nextafter(log_min - np.max(log_spectrum), -np.inf)
    high = np.nextafter(log_max - np.min(log_spectrum), np.inf)

    def evaluate(log_beta):
        bounded = bounded_log_spectrum(log_spectrum, log_beta, log_min, log_max, log_aspect)
        complexity = logsumexp(log_areas + 0.5 * bounded.sum(axis=1))
        return bounded, complexity

    for beta in (0.0, low, high):
        bounded, complexity = evaluate(beta)
        # Prefer the original scale, but return exact saturation for endpoint targets.
        endpoint = abs(lower_error) <= slack or abs(upper_error) <= slack
        if abs(complexity - log_target) <= rtol and (
            not endpoint or np.all(bounded == (log_min if abs(lower_error) <= slack else log_max))
        ):
            return bounded, float(beta), float(np.exp(complexity)), interval, 0
    for iteration in range(1, parameters.bisection_steps + 1):
        beta = low + 0.5 * (high - low)
        bounded, complexity = evaluate(beta)
        error = complexity - log_target
        if abs(error) <= rtol:
            return bounded, float(beta), float(np.exp(complexity)), interval, iteration
        if beta == low or beta == high:
            break
        if error < 0:
            low = beta
        else:
            high = beta
    raise ArithmeticError("bounded metric complexity did not reach the requested tolerance")


def construct_metric(mesh: Mesh, hessian, parameters: MetricParameters):
    """Construct K, M_star and M_b from nodal Hessians (N,2,2)."""
    eigenvalues, basis, curvature = regularise_curvature(hessian, parameters.epsilon)
    if len(eigenvalues) != len(mesh.coordinates):
        raise ValueError("expected one Hessian per mesh vertex")
    log_optimal, log_s = optimal_log_spectrum(curvature, mesh.nodal_areas, parameters.p,
                                               parameters.complexity)
    log_bounded, log_beta, c_b, interval, iterations = recover_complexity(
        log_optimal, mesh.nodal_areas, parameters
    )
    with np.errstate(over="ignore", under="ignore"):
        optimal, bounded = np.exp(log_optimal), np.exp(log_bounded)
    tensors = []
    for spectrum in (curvature, optimal, bounded):
        if np.any(spectrum <= 0) or not np.isfinite(spectrum).all():
            raise ArithmeticError("metric or curvature spectrum exceeds floating-point range")
        matrix = tensor(basis, spectrum)
        if not np.isfinite(matrix).all() or np.any(np.linalg.eigvalsh(matrix) <= 0):
            raise ArithmeticError("metric or curvature tensor cannot retain SPD in Cartesian coordinates")
        tensors.append(matrix)
    c_opt = np.exp(logsumexp(np.log(mesh.nodal_areas) + 0.5 * log_optimal.sum(axis=1)))
    return MetricResult(eigenvalues, basis, curvature, *tensors, optimal, bounded,
                        float(c_opt), c_b, interval, log_s, log_beta, iterations)
