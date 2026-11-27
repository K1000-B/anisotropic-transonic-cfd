"""Public entry point for a single solver-independent adaptation."""

from dataclasses import dataclass
from time import perf_counter

import numpy as np

from .field import P1MetricField
from .mesh import Mesh, readonly
from .metric import MetricParameters, MetricResult, Recovery, construct_metric, recover_hessian
from .remesh import GeometryConstraints, RemeshDiagnostics, RemeshParameters, remesh


@dataclass(frozen=True)
class AdaptationResult:
    mesh: Mesh
    density: np.ndarray
    recovery: Recovery
    metric: MetricResult
    diagnostics: RemeshDiagnostics
    source_vertex_indices: np.ndarray
    timings: dict[str, float]


def adapt(
    mesh: Mesh,
    density,
    metric: MetricParameters,
    remeshing=RemeshParameters(),
    constraints=GeometryConstraints(),
):
    """Adapt (source Mesh, aligned nodal density of shape (N,)) in one cycle.

    The returned mesh has fresh contiguous output indices. The result retains
    source-vertex indices and fields for diagnostics. Density stays on the
    source grid and is available for the caller's chosen transfer procedure.
    """
    start = perf_counter()
    density = readonly(np.asarray(density, dtype=float))
    recovery = recover_hessian(mesh, density)
    recovered = perf_counter()
    metrics = construct_metric(mesh, recovery.hessian, metric)
    field = P1MetricField(mesh, metrics.bounded, tolerance=remeshing.location_tolerance)
    constructed = perf_counter()
    adapted, diagnostics, origin = remesh(field, remeshing, constraints)
    return AdaptationResult(
        adapted, density, recovery, metrics, diagnostics, origin,
        {"hessian_seconds": recovered - start, "metric_and_index_seconds": constructed - recovered,
         "remesh_seconds": perf_counter() - constructed},
    )
