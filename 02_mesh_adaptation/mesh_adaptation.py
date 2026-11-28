
"""Hessian-based anisotropic adaptation of a 2D triangular SU2 mesh.

Input: NPZ produced by extract_solution.py (seven arrays).
Output: a new SU2 ASCII mesh, with its original boundary markers.

The source mesh supports a fixed, continuous P1 metric field. The working
mesh is modified independently by edge splits, collapses, diagonal swaps
and vertex relocation. No CFD solution is interpolated or recomputed here.

Only NumPy and the Python standard library are required.
"""

import argparse
from collections import defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np


# Numerical safeguards, not physical metric parameters.
COMPLEXITY_RTOL = 1e-10
BARY_TOL = 1e-10
QUALITY_GAIN = 1e-6
LOCAL_QUALITY = 0.08
FINAL_QUALITY = 0.20
MIDPOINT_RTOL = 1e-7
PROTECTED_TAGS = ("airfoil", "farfield")


# ============================================================
# 1. INPUT AND SOURCE MESH GEOMETRY
# ============================================================

def load_input(input_path):
    """Read the seven aligned NumPy arrays without changing their numbering."""

    with np.load(input_path, allow_pickle=False) as data:
        points = data['points'].astype(float)
        triangles = data['triangles']
        scalar = data['scalar'].astype(float)
        boundary_edges = data['boundary_edges']
        boundary_tags = data['boundary_tags'].astype(str)
        boundary_nodes = data['boundary_nodes']
        field_name = data['field_name'].item()

    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 3:
        raise ValueError('Expected a nonempty (N, 2) points array.')

    if scalar.shape != (len(points),):
        raise ValueError('Invalid nodal scalar array.')

    if not np.isfinite(points).all() or not np.isfinite(scalar).all():
        raise ValueError('Non-finite coordinates or scalar values.')

    for array, width in ((triangles, 3), (boundary_edges, 2)):
        if (array.ndim != 2 or array.shape[1] != width or not len(array)
                or not np.issubdtype(array.dtype, np.integer)
                or np.any(array < 0) or np.any(array >= len(points))):
            raise ValueError('Invalid mesh connectivity.')

    if boundary_tags.shape != (len(boundary_edges),):
        raise ValueError('Boundary tags and edges are not aligned.')

    if np.any(boundary_tags == ''):
        raise ValueError('Empty boundary marker.')

    if (boundary_nodes.ndim != 1
            or not np.issubdtype(boundary_nodes.dtype, np.integer)
            or np.any(boundary_nodes < 0)
            or np.any(boundary_nodes >= len(points))):
        raise ValueError('Invalid boundary_nodes array.')

    # The extractor must protect airfoil and farfield, but not outlet-only
    # vertices. Reject a stale NPZ rather than silently changing this policy.

    expected = np.unique(
        boundary_edges[np.isin(boundary_tags, PROTECTED_TAGS)]
    )

    if not np.array_equal(np.unique(boundary_nodes), expected):
        raise ValueError(
            'Regenerate the NPZ: boundary_nodes must contain only '
            'the unique airfoil and farfield vertex indices.'
        )

    return (
        points,
        triangles.copy(),
        scalar,
        boundary_edges,
        boundary_tags,
        expected,
        field_name
    )


def triangle_vertices(T, triangles):
    """Return the three vertex indices i, j, l of source triangle T."""

    return triangles[T]


def double_area(a, b, c):
    """Twice the signed physical area; positive for CCW triangles."""

    u, v = b - a, c - a

    return float(u[0] * v[1] - u[1] * v[0])


def triangle_geometry(points, triangles):
    """Build B_T and |T| for every source triangle, exactly as in the theory."""

    x = points[triangles]

    B = np.stack(
        (x[:, 1] - x[:, 0], x[:, 2] - x[:, 0]),
        axis=1
    )

    signed = np.linalg.det(B)

    if not np.isfinite(signed).all() or np.any(signed == 0):
        raise ValueError('Degenerate source triangle.')

    # Orient the working/source connectivity consistently,
    # without moving nodes.

    reverse = signed < 0

    triangles[reverse] = triangles[reverse][:, [0, 2, 1]]

    x = points[triangles]

    B = np.stack(
        (x[:, 1] - x[:, 0], x[:, 2] - x[:, 0]),
        axis=1
    )

    return B, 0.5 * np.abs(signed)


def incident_triangles(triangles, n_vertices):
    """Build vertex-to-triangle adjacency once for differential recovery."""

    incident = [[] for _ in range(n_vertices)]

    for T, (i, j, l) in enumerate(triangles):

        for vertex in (i, j, l):

            incident[vertex].append(T)

    return [
        np.asarray(ids, dtype=int)
        for ids in incident
    ]


# ============================================================
# 2. DOUBLE P1 DIFFERENTIAL RECONSTRUCTION
# ============================================================

def triangle_gradient(T, triangles, scalar, B):
    """Solve B_T g_T = [phi_j - phi_i, phi_l - phi_i] for one triangle."""

    i, j, l = triangle_vertices(T, triangles)

    return np.linalg.solve(
        B[T],
        [scalar[j] - scalar[i], scalar[l] - scalar[i]]
    )


def element_gradients(triangles, scalar, B):
    """The same linear solve for all triangles, using NumPy batch operations."""

    # Retrieve the scalar values of the three vertices
    # of every triangle.

    f = scalar[triangles]

    # Construct the edge differences of the scalar field.

    rhs = f[:, 1:] - f[:, :1]

    # Solve all independent 2x2 systems simultaneously.

    return np.linalg.solve(B, rhs[..., None])[..., 0]


def recover_at_nodes(element_values, incident, areas, area_sums):
    """Area-weighted nodal average of an elementwise vector field."""

    result = np.empty(
        (len(incident), element_values.shape[1])
    )

    for i, ids in enumerate(incident):

        result[i] = np.sum(
            areas[ids, None] * element_values[ids],
            axis=0
        ) / area_sums[i]

    return result


def recover_hessian(
    triangles,
    scalar,
    B,
    incident,
    areas,
    area_sums
):
    """Recover gradients, differentiate each component, then symmetrise D."""

    # First P1 reconstruction: gradient of the scalar field.

    g_T = element_gradients(
        triangles,
        scalar,
        B
    )

    # Recover elemental gradients at the mesh vertices.

    g_i = recover_at_nodes(
        g_T,
        incident,
        areas,
        area_sums
    )

    # Second P1 reconstruction: gradient of each component
    # of the recovered nodal gradient.
    #
    # Row 0 of D is grad(g_x).
    # Row 1 of D is grad(g_y).

    rows = []

    for component in range(2):

        dg_T = element_gradients(
            triangles,
            g_i[:, component],
            B
        )

        rows.append(
            recover_at_nodes(
                dg_T,
                incident,
                areas,
                area_sums
            )
        )

    D = np.stack(rows, axis=1)

    # Symmetrise the Hessian to average the two estimates
    # of the mixed derivatives.

    return 0.5 * (
        D + D.transpose(0, 2, 1)
    )


# ============================================================
# 3. CURVATURE AND ANALYTICAL OPTIMAL METRIC
# ============================================================

def regularised_curvature(H, epsilon):
    """
    K = Q diag(max(abs(lambda), epsilon)) Q^T
    det(K) = k1*k2
    """

    # Spectral decomposition of the symmetric Hessian.

    eigenvalues, Q = np.linalg.eigh(H)

    # Take the absolute value of the Hessian eigenvalues
    # and impose the strictly positive curvature floor.

    k = np.maximum(
        np.abs(eigenvalues),
        epsilon
    )

    # Reconstruct the regularised curvature tensors.

    K = (
        Q * k[:, None, :]
    ) @ Q.transpose(0, 2, 1)

    # Compute the determinant directly from the eigenvalues.

    det_K = k[:, 0] * k[:, 1]

    if (
        not np.isfinite(K).all()
        or not np.isfinite(det_K).all()
        or np.any(det_K <= 0)
    ):
        raise ValueError(
            'Invalid curvature spectrum; check epsilon and field scale.'
        )

    return K, det_K


def optimal_metric(A, K, det_K, target_nodes, p):
    """Apply the analytical result, using complexity as the actual budget."""

    # Convert the indicative vertex budget into the target
    # metric complexity.

    C_star = (
        np.sqrt(3.0) / 2.0
    ) * target_nodes

    # Compute the global scaling factor:
    #
    # s_star = C_star /
    #          sum_i(A_i * det(K_i)^(p/(2*(p+1))))

    s_star = C_star / np.sum(
        A * det_K ** (p / (2.0 * (p + 1)))
    )

    # Construct the optimal metric:
    #
    # M_star_i =
    #     s_star * det(K_i)^(-1/(2*(p+1))) * K_i

    M_star = (
        s_star * det_K ** (-1.0 / (2.0 * (p + 1)))
    )[:, None, None] * K

    if (
        not np.isfinite(M_star).all()
        or not np.isfinite(s_star)
        or s_star <= 0
    ):
        raise ValueError('Invalid optimal metric scaling.')

    return M_star, C_star, s_star


def complexity(A, M):
    """Discrete complexity: sum_i A_i sqrt(det(M_i))."""

    eig = np.linalg.eigvalsh(M)

    if np.any(eig <= 0) or not np.isfinite(eig).all():
        raise ValueError('Metric is not positive definite.')

    return float(
        np.sum(
            A
            * np.sqrt(eig[:, 0])
            * np.sqrt(eig[:, 1])
        )
    )


def bounded_metric(
    M_star,
    A,
    C_star,
    h_min,
    h_max,
    max_aspect_ratio
):
    """Apply spectral bounds and restore complexity by log-space bisection."""

    # The target element sizes and metric eigenvalues satisfy:
    #
    # h = 1 / sqrt(lambda)
    #
    # Consequently:
    #
    # 1 / h_max^2 <= lambda <= 1 / h_min^2

    low = 1.0 / h_max**2
    high = 1.0 / h_min**2

    area = np.sum(A)

    # Verify the feasibility of the prescribed complexity.

    C_min = area * low
    C_max = area * high

    if (
        C_star < C_min * (1 - COMPLEXITY_RTOL)
        or C_star > C_max * (1 + COMPLEXITY_RTOL)
    ):
        raise ValueError(
            f'Infeasible complexity {C_star:g}; '
            f'interval [{C_min:g}, {C_max:g}].'
        )

    # Diagonalise the optimal metric once.
    #
    # Positive scalar multiplication modifies the eigenvalues
    # but does not modify the eigenvectors.

    eigenvalues, Q = np.linalg.eigh(M_star)

    if np.any(eigenvalues <= 0) or not np.isfinite(eigenvalues).all():
        raise ValueError('Invalid optimal metric spectrum.')

    def evaluate(log_beta):

        # Scale BEFORE bounding.
        # The directions of the eigenvectors remain fixed.

        values = np.exp(
            np.clip(
                log_beta + np.log(eigenvalues),
                np.log(low),
                np.log(high)
            )
        )

        # The eigenvalues are sorted in ascending order.
        #
        # lambda_minus = values[:, 0]
        # lambda_plus  = values[:, 1]
        #
        # Enforce:
        #
        # sqrt(lambda_plus / lambda_minus) <= A_max

        values[:, 0] = np.maximum(
            values[:, 0],
            values[:, 1] / max_aspect_ratio**2
        )

        # Compute the complexity of the bounded metric.

        C = np.sum(
            A
            * np.sqrt(values[:, 0])
            * np.sqrt(values[:, 1])
        )

        return C, values

    # Construct the initial bisection interval.
    #
    # At the lower endpoint, all metric eigenvalues
    # reach their minimum allowed value.
    #
    # At the upper endpoint, all metric eigenvalues
    # reach their maximum allowed value.

    left = np.log(low) - np.log(np.max(eigenvalues))

    right = np.log(high) - np.log(np.min(eigenvalues))

    # Endpoints also handle an isotropic fixed-size metric
    # and saturated budgets.

    for trial in (left, right):

        C, values = evaluate(trial)

        if abs(C - C_star) / C_star <= COMPLEXITY_RTOL:

            return (
                (Q * values[:, None, :]) @ Q.transpose(0, 2, 1),
                np.exp(trial)
            )

    # Restore the complexity using bisection in log(beta).

    for _ in range(100):

        trial = 0.5 * (left + right)

        C, values = evaluate(trial)

        if abs(C - C_star) / C_star <= COMPLEXITY_RTOL:

            return (
                (Q * values[:, None, :]) @ Q.transpose(0, 2, 1),
                np.exp(trial)
            )

        if C < C_star:

            left = trial

        else:

            right = trial

    raise RuntimeError(
        'Bounded metric complexity bisection failed.'
    )


# ============================================================
# 4. FIXED PIECEWISE-P1 METRIC FIELD
# ============================================================

def build_metric_field(points, triangles, M):
    """
    Return metric evaluation, exact edge length and metric
    midpoint functions.

    A simple uniform grid indexes source triangle bounding boxes.

    Segment intersections are computed in barycentric coordinates.

    The source mesh is never changed, including when a working
    vertex is relocated.
    """

    # --------------------------------------------------------
    # 4.1. Prepare barycentric transformations
    # --------------------------------------------------------

    x = points[triangles]

    origins = x[:, 0]

    inverse = np.linalg.inv(
        np.stack(
            (
                x[:, 1] - origins,
                x[:, 2] - origins
            ),
            axis=2
        )
    )

    # Store the metric tensors of the three vertices
    # of every source triangle.

    nodal_metrics = M[triangles]

    # --------------------------------------------------------
    # 4.2. Construct a spatial search grid
    # --------------------------------------------------------

    # A spatial grid avoids searching through every source
    # triangle for each metric evaluation.

    lower = points.min(axis=0)

    extent = points.max(axis=0) - lower

    n = max(
        1,
        int(np.sqrt(len(triangles) / 8))
    )

    step = extent / n

    cells = defaultdict(list)

    def cell(x):

        return np.clip(
            np.floor((x - lower) / step).astype(int),
            0,
            n - 1
        )

    # Register each source triangle in all cells intersected
    # by its axis-aligned bounding box.

    for T, vertices in enumerate(x):

        first = cell(vertices.min(axis=0))

        last = cell(vertices.max(axis=0))

        for i in range(first[0], last[0] + 1):

            for j in range(first[1], last[1] + 1):

                cells[i, j].append(T)

    # --------------------------------------------------------
    # 4.3. Retrieve candidate source triangles
    # --------------------------------------------------------

    def candidates(a, b):

        first = cell(np.minimum(a, b))

        last = cell(np.maximum(a, b))

        ids = set()

        for i in range(first[0], last[0] + 1):

            for j in range(first[1], last[1] + 1):

                ids.update(
                    cells.get((i, j), ())
                )

        return np.asarray(
            sorted(ids),
            dtype=int
        )

    # --------------------------------------------------------
    # 4.4. Compute barycentric coordinates
    # --------------------------------------------------------

    def barycentric(a, ids):

        z = np.einsum(
            'nij,nj->ni',
            inverse[ids],
            a - origins[ids]
        )

        return np.column_stack(
            (
                1.0 - z.sum(axis=1),
                z
            )
        )

    # --------------------------------------------------------
    # 4.5. Evaluate the continuous metric at a physical point
    # --------------------------------------------------------

    @lru_cache(maxsize=100000)
    def metric_cached(a):

        a = np.asarray(a)

        ids = candidates(a, a)

        w = barycentric(a, ids)

        # A point belongs to a triangle when its three
        # barycentric coordinates are nonnegative.

        inside = np.flatnonzero(
            np.all(
                w >= -BARY_TOL,
                axis=1
            )
        )

        if not len(inside):
            raise ValueError('Point outside the source mesh.')

        k = inside[0]

        # Correct negligible numerical violations of
        # barycentric nonnegativity.

        weights = np.maximum(
            w[k],
            0.0
        )

        weights /= weights.sum()

        # P1 interpolation:
        #
        # M(x) = N_i*M_i + N_j*M_j + N_l*M_l

        return np.einsum(
            'i,ijk->jk',
            weights,
            nodal_metrics[ids[k]]
        )

    def metric_at(a):

        return metric_cached(tuple(a))

    # --------------------------------------------------------
    # 4.6. Exact integration of the metric length
    # --------------------------------------------------------

    def integrate(q0, q1, dt):
        """
        Integrate sqrt(q(t)) exactly over an interval where
        q(t) is affine.

        q(t) = d^T M(x(t)) d

        Integral:

            2*dt/3 * (q0 + sqrt(q0*q1) + q1)
                   / (sqrt(q0) + sqrt(q1))
        """

        if min(q0, q1) <= 0 or not np.isfinite(q0 + q1):
            raise ValueError('Invalid edge metric integrand.')

        a, b = np.sqrt(q0), np.sqrt(q1)

        return (
            (2.0 * dt / 3.0)
            * (q0 + a*b + q1)
            / (a+b)
        )

    # --------------------------------------------------------
    # 4.7. Decompose an edge into source triangle intervals
    # --------------------------------------------------------

    @lru_cache(maxsize=100000)
    def edge_pieces(a, b):
        """
        Find non-overlapping intervals covering the segment
        in the source mesh.

        Inside each triangle, each barycentric weight is affine
        in t. Its nonnegativity defines a linear constraint on t.

        Intersect the three constraints with [0,1].

        Shared edges may give duplicate intervals; these are
        covered only once.

        No metric extrapolation is permitted.
        """

        a, b = np.asarray(a), np.asarray(b)

        d = b - a

        if not np.any(d):
            raise ValueError('Zero-length edge.')

        ids = candidates(a, b)

        # Compute barycentric coordinates at both endpoints.

        w0, w1 = barycentric(a, ids), barycentric(b, ids)

        dw = w1 - w0

        # Every source triangle initially has the complete
        # parameter interval [0,1].

        lo = np.zeros(len(ids))

        hi = np.ones(len(ids))

        valid = np.ones(
            len(ids),
            dtype=bool
        )

        # Intersect the interval with each barycentric
        # nonnegativity constraint.

        for k in range(3):

            w, delta = w0[:, k], dw[:, k]

            up = delta > BARY_TOL

            down = delta < -BARY_TOL

            flat = ~(up | down)

            valid &= ~(flat & (w < -BARY_TOL))

            lo[up] = np.maximum(
                lo[up],
                -w[up] / delta[up]
            )

            hi[down] = np.minimum(
                hi[down],
                -w[down] / delta[down]
            )

        valid &= hi > lo

        # Nodal q values are positive.
        # q varies affinely in each source triangle.

        qnodes = np.einsum(
            'i,nkij,j->nk',
            d,
            nodal_metrics[ids],
            d
        )

        q0 = np.sum(w0*qnodes, axis=1)

        dq = np.sum(dw*qnodes, axis=1)

        # Sort the intervals according to their position
        # along the edge.

        order = sorted(
            np.flatnonzero(valid),
            key=lambda k: (lo[k], -hi[k])
        )

        pieces = []

        covered = 0.0

        for k in order:

            end = min(1.0, hi[k])

            # Avoid integrating the same interval twice.

            if end <= covered:
                continue

            if lo[k] > covered + BARY_TOL:
                raise ValueError('Edge leaves the source domain.')

            start = covered

            # Evaluate the affine integrand at both ends
            # of the current interval.

            qa = q0[k] + start*dq[k]

            qb = q0[k] + end*dq[k]

            L = integrate(
                qa,
                qb,
                end-start
            )

            pieces.append(
                (start, end, qa, qb, L)
            )

            covered = end

            if covered >= 1.0:
                break

        if covered < 1.0 - BARY_TOL or not pieces:
            raise ValueError(
                'Edge not covered by the source mesh.'
            )

        return tuple(pieces)

    # --------------------------------------------------------
    # 4.8. Compute the complete metric length of an edge
    # --------------------------------------------------------

    def metric_length(a, b):

        # Canonical endpoint order improves cache reuse
        # and ensures numerical symmetry.

        a, b = tuple(a), tuple(b)

        if a == b:
            return 0.0

        return sum(
            piece[4]
            for piece in edge_pieces(
                *sorted((a, b))
            )
        )

    # --------------------------------------------------------
    # 4.9. Compute the metric midpoint
    # --------------------------------------------------------

    def metric_midpoint(a, b):
        """
        Find the point dividing the integrated metric length
        into two equal halves.

        The returned point P satisfies:

            L_M(A, P) = L_M(P, B)

        It does not necessarily correspond to the Euclidean
        midpoint of the edge.
        """

        a, b = np.asarray(a), np.asarray(b)

        pieces = edge_pieces(
            tuple(a),
            tuple(b)
        )

        target = 0.5 * sum(
            piece[4]
            for piece in pieces
        )

        remaining = target

        for start, end, qa, qb, L in pieces:

            if remaining > L:

                remaining -= L

                continue

            # Invert the cumulative integral on this one
            # affine interval using bisection.

            left, right = 0.0, 1.0

            for _ in range(30):

                fraction = 0.5 * (left + right)

                qmid = qa + fraction * (qb - qa)

                value = integrate(
                    qa,
                    qmid,
                    fraction * (end-start)
                )

                if abs(value - remaining) <= MIDPOINT_RTOL * target:
                    break

                if value < remaining:

                    left = fraction

                else:

                    right = fraction

            t = start + fraction * (end-start)

            return a + t * (b-a)

        raise RuntimeError('Metric midpoint not found.')

    return metric_at, metric_length, metric_midpoint


# ============================================================
# 5. WORKING MESH AND LOCAL OPERATIONS
# ============================================================

def edge_key(a, b):
    """Return a canonical unoriented edge identifier."""

    return (min(a, b), max(a, b))


def triangle_edges(t):
    """Return the three unoriented edges of a triangle."""

    i, j, l = t

    return (
        edge_key(i, j),
        edge_key(j, l),
        edge_key(l, i)
    )


def metric_quality(coords, metric_at):
    """
    Mean ratio using the average of the three vertex
    metric tensors.

    Q_M = 4*sqrt(3)*A_M / sum(L_M^2)

    A value of one corresponds to an equilateral triangle
    in metric space.
    """

    a, b, c = coords

    area2 = double_area(a, b, c)

    if area2 <= 0:
        return 0.0

    # Approximate the metric inside the triangle by the
    # arithmetic mean of its three vertex metrics.

    M = (
        metric_at(a)
        + metric_at(b)
        + metric_at(c)
    ) / 3.0

    # Compute the three physical edge vectors.

    d = np.array([
        b-a,
        c-b,
        a-c
    ])

    # Sum of squared edge lengths in the averaged metric.

    den = np.einsum(
        'ni,ij,nj->',
        d,
        M,
        d
    )

    return float(
        2
        * np.sqrt(3)
        * area2
        * np.sqrt(np.linalg.det(M))
        / den
    )


def adapt_mesh(
    source_points,
    source_triangles,
    boundary_edges,
    boundary_tags,
    boundary_nodes,
    metric_at,
    metric_length,
    metric_midpoint,
    collapse_threshold,
    split_threshold,
    max_iterations,
    diagnostics=None
):
    """
    Perform successive local edits using a fixed source metric.

    Dictionaries provide stable integer IDs while vertices and
    triangles are added or removed.

    The adjacency tables are updated only on the edited patch;
    they are not rebuilt for every individual operation.
    """

    # --------------------------------------------------------
    # 5.1. Initialise the working mesh
    # --------------------------------------------------------

    # Copy the source coordinates to ensure that the
    # original metric interpolation support remains unchanged.

    points = {
        i: x.copy()
        for i, x in enumerate(source_points)
    }

    triangles = {}

    edge_owners = defaultdict(set)

    vertex_owners = defaultdict(set)

    # Store oriented boundary records, indexed by an
    # unoriented edge key.

    boundaries = {}

    for (a, b), tag in zip(boundary_edges, boundary_tags):

        key = edge_key(int(a), int(b))

        if key in boundaries:
            raise ValueError('Duplicate marked boundary edge.')

        boundaries[key] = (
            int(a),
            int(b),
            str(tag)
        )

    original_boundaries = boundaries.copy()

    protected = set(map(int, boundary_nodes))

    # Unrecognised markers are frozen conservatively.
    # Only outlet is editable.

    for a, b, tag in boundaries.values():

        if tag != 'outlet':

            protected.update((a, b))

    # Store boundary adjacency independently from the
    # complete triangle adjacency.

    boundary_neighbors = defaultdict(set)

    for a, b, _ in boundaries.values():

        boundary_neighbors[a].add(b)

        boundary_neighbors[b].add(a)

    next_vertex = len(points)

    next_triangle = 0

    # --------------------------------------------------------
    # 5.2. Triangle connectivity management
    # --------------------------------------------------------

    def add_triangle(t):

        nonlocal next_triangle

        T = next_triangle

        next_triangle += 1

        triangles[T] = tuple(t)

        # Associate the triangle with its three vertices.

        for i in t:

            vertex_owners[i].add(T)

        # Associate the triangle with its three edges.

        for e in triangle_edges(t):

            edge_owners[e].add(T)

    def remove_triangle(T):

        t = triangles.pop(T)

        for i in t:

            vertex_owners[i].remove(T)

        for e in triangle_edges(t):

            edge_owners[e].remove(T)

            if not edge_owners[e]:

                del edge_owners[e]

    # Initialise triangle connectivity and adjacency.

    for t in source_triangles:

        add_triangle(
            tuple(map(int, t))
        )

    # --------------------------------------------------------
    # 5.3. Boundary and topology checks
    # --------------------------------------------------------

    def check_boundary_topology():

        if any(
            len(owners) not in (1, 2)
            for owners in edge_owners.values()
        ):
            raise ValueError('Non-manifold edge.')

        if {
            e for e, owners in edge_owners.items()
            if len(owners) == 1
        } != set(boundaries):

            raise ValueError(
                'Boundary markers do not match mesh connectivity.'
            )

    check_boundary_topology()

    if len({
        tuple(sorted(t))
        for t in triangles.values()
    }) != len(triangles):

        raise ValueError('Duplicate source triangles.')

    # --------------------------------------------------------
    # 5.4. Local geometric utilities
    # --------------------------------------------------------

    def neighbors(i):
        """Return all vertices connected to vertex i."""

        return {
            j
            for T in vertex_owners[i]
            for j in triangles[T]
            if j != i
        }

    def L(e):
        """Return the metric length of a working edge."""

        return metric_length(
            points[e[0]],
            points[e[1]]
        )

    def coordinates(t, changes):
        """
        Retrieve the coordinates of the vertices of t.

        Candidate coordinates take precedence over the
        current coordinates of the working mesh.
        """

        return np.asarray([
            changes[i] if i in changes else points[i]
            for i in t
        ])

    def Q(t, changes=None):

        return metric_quality(
            coordinates(t, changes or {}),
            metric_at
        )

    def orient_new(t, changes):

        # New split/swap triangles can be listed in either
        # order.
        #
        # Collapse and relocation must instead preserve
        # the old orientation strictly.

        a, b, c = coordinates(t, changes)

        return (
            t
            if double_area(a, b, c) > 0
            else (t[0], t[2], t[1])
        )

    # --------------------------------------------------------
    # 5.5. Validate candidate local mesh modifications
    # --------------------------------------------------------

    def patch_valid(
        old_ids,
        new_triangles,
        changes,
        improve=False,
        guard_lengths=False
    ):
        """
        Check a candidate local replacement before committing
        any change.

        The checks cover:

            - Duplicate triangles.
            - Invalid vertex connectivity.
            - Non-manifold edges.
            - Degenerate or inverted triangles.
            - Conservation of local physical area.
            - Metric triangle quality.
            - Excessively long metric edges.
        """

        old_ids = set(old_ids)

        if not new_triangles:
            return False, 0.0

        old = [
            triangles[T]
            for T in old_ids
        ]

        keys = [
            tuple(sorted(t))
            for t in new_triangles
        ]

        if len(set(keys)) != len(keys):
            return False, 0.0

        # Ensure edge incidence after replacement remains
        # manifold and prevent duplicates with triangles
        # outside the affected patch.

        counts = defaultdict(int)

        for t in old:

            for e in triangle_edges(t):

                counts[e] -= 1

        for t in new_triangles:

            if len(set(t)) != 3:
                return False, 0.0

            for T in vertex_owners.get(t[0], set()) - old_ids:

                if tuple(sorted(triangles[T])) == tuple(sorted(t)):
                    return False, 0.0

            for e in triangle_edges(t):

                counts[e] += 1

        if any(
            not 0 <= len(edge_owners.get(e, ())) + delta <= 2
            for e, delta in counts.items()
        ):
            return False, 0.0

        old_area = sum(
            double_area(*coordinates(t, {}))
            for t in old
        )

        new_area = 0.0

        try:

            # Check the orientation and degeneracy of every
            # candidate triangle.

            for t in new_triangles:

                x = coordinates(t, changes)

                area = double_area(*x)

                scale = max(
                    np.dot(x[j]-x[i], x[j]-x[i])
                    for i, j in ((0, 1), (1, 2), (2, 0))
                )

                if area <= 64 * np.finfo(float).eps * scale:
                    return False, 0.0

                new_area += area

            # An unchanged patch boundary must enclose
            # the same physical area.

            if abs(new_area - old_area) > 1e-9 * old_area:
                return False, 0.0

            # Compare the minimum metric quality before
            # and after the candidate operation.

            before = min(
                Q(t)
                for t in old
            )

            after = min(
                Q(t, changes)
                for t in new_triangles
            )

            if after < min(LOCAL_QUALITY, before) - QUALITY_GAIN:
                return False, 0.0

            # Swap and relocation operations must provide
            # a strictly positive quality improvement.

            if improve and after <= before + QUALITY_GAIN:
                return False, 0.0

            if guard_lengths:

                old_edges = {
                    e
                    for t in old
                    for e in triangle_edges(t)
                }

                new_edges = {
                    e
                    for t in new_triangles
                    for e in triangle_edges(t)
                }

                max_before = max(
                    L(e)
                    for e in old_edges
                )

                # Check newly created or geometrically
                # modified edges.

                for e in new_edges:

                    if (
                        e not in old_edges
                        or any(i in changes for i in e)
                    ):

                        a, b = coordinates(e, changes)

                        if metric_length(a, b) > (
                            max(split_threshold, max_before)
                            * (1 + QUALITY_GAIN)
                        ):
                            return False, 0.0

            # Positive oriented triangles, link checks and
            # unchanged patch boundary prevent foldovers
            # for a conforming initial mesh.

            return True, after

        except (ValueError, FloatingPointError):

            return False, 0.0

    # --------------------------------------------------------
    # 5.6. Commit a valid mesh modification
    # --------------------------------------------------------

    def commit(old_ids, new_triangles, changes):

        for T in sorted(old_ids):

            remove_triangle(T)

        points.update(changes)

        for t in new_triangles:

            add_triangle(t)

    # --------------------------------------------------------
    # 5.7. Boundary connectivity management
    # --------------------------------------------------------

    def remove_boundary(e):

        a, b, _ = boundaries.pop(e)

        boundary_neighbors[a].remove(b)

        boundary_neighbors[b].remove(a)

    def add_boundary(a, b, tag):

        boundaries[edge_key(a, b)] = (
            a,
            b,
            tag
        )

        boundary_neighbors[a].add(b)

        boundary_neighbors[b].add(a)

    # --------------------------------------------------------
    # 5.8. EDGE SPLIT
    # --------------------------------------------------------

    def split(e):
        """
        Insert a new vertex at the metric midpoint of an edge.

        Interior edge:

            2 triangles -> 4 triangles

        Outlet boundary edge:

            1 triangle -> 2 triangles

        Airfoil and farfield edges cannot be subdivided.
        """

        nonlocal next_vertex

        a, b = e

        record = boundaries.get(e)

        # Protected boundary edges cannot be subdivided.

        if record is not None and record[2] != 'outlet':
            return False, None

        # Find all triangles sharing the edge.

        owners = set(edge_owners[e])

        # Assign an identifier to the candidate new vertex.

        P = next_vertex

        # Place the new vertex at the metric midpoint.

        changes = {
            P: metric_midpoint(points[a], points[b])
        }

        new = []

        # Replace every incident triangle (a,b,c)
        # by two triangles:
        #
        #     (a,P,c)
        #     (P,b,c)

        for T in sorted(owners):

            c = next(
                i for i in triangles[T]
                if i not in e
            )

            new.extend([
                orient_new((a, P, c), changes),
                orient_new((P, b, c), changes)
            ])

        # Reject the operation if the resulting patch
        # is invalid.

        ok, _ = patch_valid(
            owners,
            new,
            changes
        )

        if not ok:
            return False, None

        # Apply the subdivision.

        commit(
            owners,
            new,
            changes
        )

        next_vertex += 1

        # If the edge belongs to the outlet, preserve
        # its physical marker on both new boundary edges.

        if record is not None:

            u, v, tag = record

            remove_boundary(e)

            add_boundary(u, P, tag)

            add_boundary(P, v, tag)

        return True, P

    # --------------------------------------------------------
    # 5.9. Identify movable outlet vertices
    # --------------------------------------------------------

    def straight_outlet_vertex(i):
        """
        Return the two adjacent outlet nodes if i lies
        on a straight boundary run.

        A corner or marker junction cannot move, even if
        not explicitly protected, because its motion would
        change the physical polygon.
        """

        nbr = sorted(
            boundary_neighbors.get(i, ())
        )

        if len(nbr) != 2 or i in protected:
            return None

        a, b = nbr

        if any(
            boundaries[edge_key(i, j)][2] != 'outlet'
            for j in nbr
        ):
            return None

        u = points[a] - points[i]

        v = points[b] - points[i]

        scale = np.linalg.norm(u) * np.linalg.norm(v)

        # Both adjacent edges must be collinear.
        #
        # The negative dot product ensures that the
        # neighboring vertices lie on opposite sides of i.

        if (
            scale == 0
            or abs(u[0]*v[1] - u[1]*v[0]) > 1e-12*scale
            or np.dot(u, v) >= 0
        ):
            return None

        return a, b

    # --------------------------------------------------------
    # 5.10. EDGE COLLAPSE
    # --------------------------------------------------------

    def collapse(e, fresh):
        """
        Merge both endpoints at the metric midpoint,
        not at an old endpoint.

        The resulting vertex retains one identifier, but
        its coordinates are replaced by the metric midpoint.
        """

        a, b = e

        # Protected vertices must never be moved or removed.
        #
        # Fresh vertices are protected temporarily to avoid
        # immediately reversing a subdivision.

        if (
            a in protected
            or b in protected
            or a in fresh
            or b in fresh
        ):
            return False

        is_boundary = e in boundaries

        if is_boundary:

            # Only straight outlet edges may collapse.
            #
            # Both endpoints must be movable and their two
            # outside neighbors must remain distinct.

            na = straight_outlet_vertex(a)

            nb = straight_outlet_vertex(b)

            if na is None or nb is None:
                return False

            left = next(i for i in na if i != b)

            right = next(i for i in nb if i != a)

            if left == right or edge_key(left, right) in boundaries:
                return False

        elif boundary_neighbors.get(a) or boundary_neighbors.get(b):

            # Do not merge an interior node with a boundary
            # node: their midpoint would generally move
            # the domain boundary.

            return False

        owners = set(edge_owners[e])

        # Link condition:
        #
        # The common neighbors of a and b must be exactly
        # the vertices opposite their shared edge.

        opposite = {
            i
            for T in owners
            for i in triangles[T]
            if i not in e
        }

        if neighbors(a) & neighbors(b) != opposite:
            return False

        # Retrieve the complete patch affected by
        # the contraction.

        patch = vertex_owners[a] | vertex_owners[b]

        # Move the surviving vertex a to the metric
        # midpoint of the original edge.

        changes = {
            a: metric_midpoint(
                points[a],
                points[b]
            )
        }

        # Replace b by a in every affected triangle.
        #
        # The triangles directly sharing the contracted
        # edge disappear.

        new = [
            tuple(
                a if i == b else i
                for i in triangles[T]
            )
            for T in sorted(patch - owners)
        ]

        # Preserve orientation, connectivity, local area
        # and acceptable metric quality.

        ok, _ = patch_valid(
            patch,
            new,
            changes,
            guard_lengths=True
        )

        if not ok:
            return False

        if is_boundary:

            # Preserve the orientation and marker of both
            # outer boundary records while replacing a-b
            # by one vertex at the midpoint.

            old_records = [
                boundaries[edge_key(left, a)],
                boundaries[edge_key(b, right)]
            ]

            for key in (
                edge_key(left, a),
                e,
                edge_key(b, right)
            ):

                remove_boundary(key)

            for u, v, tag in old_records:

                add_boundary(
                    a if u == b else u,
                    a if v == b else v,
                    tag
                )

        # Apply the modified connectivity and coordinates.

        commit(
            patch,
            new,
            changes
        )

        # Remove the unused vertex and adjacency records.

        del points[b]

        vertex_owners.pop(b, None)

        boundary_neighbors.pop(b, None)

        return True

    # --------------------------------------------------------
    # 5.11. EDGE SWAP
    # --------------------------------------------------------

    def swap(e):
        """
        Replace an interior diagonal by the opposite diagonal
        if the minimum local metric quality increases.

        Original triangles:

            (A,B,C)
            (B,A,D)

        Candidate triangles:

            (C,D,A)
            (D,C,B)
        """

        # Boundary edges must not be swapped.

        if e in boundaries or len(edge_owners[e]) != 2:
            return False

        a, b = e

        owners = set(edge_owners[e])

        T, U = sorted(owners)

        # Retrieve the vertices opposite the current edge.

        c = next(
            i for i in triangles[T]
            if i not in e
        )

        d = next(
            i for i in triangles[U]
            if i not in e
        )

        # Reject invalid or already existing diagonals.

        if c == d or edge_key(c, d) in edge_owners:
            return False

        # The two old endpoints must lie on opposite sides
        # of the new diagonal.
        #
        # Together with the old valid triangles, this
        # tests the convexity of the quadrilateral.

        if (
            double_area(points[c], points[d], points[a])
            * double_area(points[c], points[d], points[b])
            >= 0
        ):
            return False

        # Construct the two candidate triangles.

        new = [
            orient_new((c, d, a), {}),
            orient_new((d, c, b), {})
        ]

        old_length = L(e)

        new_length = L(edge_key(c, d))

        def violation(length):
            """Measure the distance from the admissible length interval."""

            return max(
                collapse_threshold - length,
                length - split_threshold,
                0.0
            )

        # The new diagonal must not increase the
        # metric length violation.

        if violation(new_length) > violation(old_length) + QUALITY_GAIN:
            return False

        # Require a strictly positive local quality gain.

        ok, _ = patch_valid(
            owners,
            new,
            {},
            improve=True
        )

        if ok:

            commit(
                owners,
                new,
                {}
            )

        return ok

    # --------------------------------------------------------
    # 5.12. VERTEX RELOCATION
    # --------------------------------------------------------

    def relocate(i):
        """
        Improve local triangle quality by moving vertex i
        towards a metric-weighted target.

        Solve:

            S*x_target = rhs

        where:

            S   = sum_j M_ij

            rhs = sum_j M_ij*x_j

        The metric M_ij is evaluated at the midpoint of
        edge (i,j).
        """

        if i in protected or not vertex_owners.get(i):
            return False

        boundary_run = None

        if boundary_neighbors.get(i):

            boundary_run = straight_outlet_vertex(i)

            if boundary_run is None:
                return False

        x = points[i]

        S = np.zeros((2, 2))

        rhs = np.zeros(2)

        # Construct the metric-weighted target position.

        for j in sorted(neighbors(i)):

            Mij = metric_at(
                0.5 * (x + points[j])
            )

            S += Mij

            rhs += Mij @ points[j]

        target = np.linalg.solve(
            S,
            rhs
        )

        if boundary_run is not None:

            # Restrict outlet relocation to its straight
            # boundary segment.

            a, b = (
                points[j]
                for j in boundary_run
            )

            d = b - a

            t = np.dot(target-a, d) / np.dot(d, d)

            target = a + np.clip(
                t,
                1e-6,
                1-1e-6
            ) * d

        owners = set(vertex_owners[i])

        # Relocation modifies the coordinates only.
        # The triangle connectivity remains unchanged.

        new = [
            triangles[T]
            for T in sorted(owners)
        ]

        # Backtracking:
        #
        # x_candidate = x + alpha*(target-x)
        #
        # If a full displacement is invalid, progressively
        # reduce alpha.

        for k in range(6):

            changes = {
                i: x + (target-x)*0.5**k
            }

            ok, _ = patch_valid(
                owners,
                new,
                changes,
                improve=True,
                guard_lengths=True
            )

            if ok:

                # Connectivity is unchanged;
                # update the coordinate only.

                points.update(changes)

                return True

        return False

    # --------------------------------------------------------
    # 5.13. Mesh statistics
    # --------------------------------------------------------

    def statistics():

        # Evaluate the metric length of every working edge.

        lengths = np.asarray([
            L(e)
            for e in edge_owners
        ])

        # Evaluate the metric quality of every working triangle.

        qualities = np.asarray([
            Q(t)
            for t in triangles.values()
        ])

        return (
            int(np.count_nonzero(lengths < collapse_threshold)),
            int(np.count_nonzero(lengths > split_threshold)),
            float(qualities.min()),
            float(np.median(qualities)),
            float(np.sqrt(np.mean((lengths-1)**2)))
        )

    # ========================================================
    # 6. MAIN REMESHING LOOP
    # ========================================================

    status = 'iteration_limit'

    # Optional observation only: never feed operation counts back into edits.
    totals = dict(split_count=0, collapse_count=0, swap_count=0,
                  relocation_count=0)

    short, long, qmin, qmed, error = statistics()

    print(
        f'Initial: nodes={len(points)}, '
        f'short={short}, long={long}, '
        f'Qmedian={qmed:.4f}',
        flush=True
    )

    for iteration in range(1, max_iterations+1):

        # Vertices inserted during the current pass.
        #
        # They must not be immediately removed by the
        # contraction phase.

        fresh = set()

        n_split = 0

        n_collapse = 0

        n_swap = 0

        n_move = 0

        # ----------------------------------------------------
        # 6.1. Split excessively long edges
        # ----------------------------------------------------

        # Sort a snapshot; recheck the edge and its length
        # before each edit.

        for e in sorted(
            edge_owners,
            key=L,
            reverse=True
        ):

            if (
                e in edge_owners
                and L(e) > split_threshold
            ):

                ok, P = split(e)

                if ok:

                    fresh.add(P)

                    n_split += 1

        # ----------------------------------------------------
        # 6.2. Collapse excessively short edges
        # ----------------------------------------------------

        for e in sorted(
            edge_owners,
            key=L
        ):

            if (
                e in edge_owners
                and L(e) < collapse_threshold
            ):

                n_collapse += int(
                    collapse(e, fresh)
                )

        # ----------------------------------------------------
        # 6.3. Improve triangle connectivity
        # ----------------------------------------------------

        for e in sorted(edge_owners):

            if e in edge_owners:

                n_swap += int(
                    swap(e)
                )

        # ----------------------------------------------------
        # 6.4. Improve vertex positions
        # ----------------------------------------------------

        for i in sorted(points):

            n_move += int(
                relocate(i)
            )

        # ----------------------------------------------------
        # 6.5. Evaluate the mesh after the current pass
        # ----------------------------------------------------

        check_boundary_topology()

        short, long, qmin, qmed, error = statistics()

        totals['split_count'] += n_split
        totals['collapse_count'] += n_collapse
        totals['swap_count'] += n_swap
        totals['relocation_count'] += n_move

        print(
            f'Pass {iteration}: '
            f'nodes={len(points)}, '
            f'split={n_split}, '
            f'collapse={n_collapse}, '
            f'swap={n_swap}, '
            f'move={n_move}, '
            f'short={short}, '
            f'long={long}, '
            f'Qmin={qmin:.4f}, '
            f'Qmedian={qmed:.4f}, '
            f'edge_RMS={error:.4f}',
            flush=True
        )

        # ----------------------------------------------------
        # 6.6. Convergence criteria
        # ----------------------------------------------------

        # The target node count is not a convergence criterion.
        # The metric complexity has already prescribed the
        # desired global resolution.

        if (
            short == 0
            and long == 0
            and qmin >= FINAL_QUALITY
        ):

            status = 'converged'

            break

        # If no operation can be applied, the mesh is
        # stationary, but not necessarily converged.

        if n_split + n_collapse + n_swap + n_move == 0:

            status = 'stationary'

            break

    # ========================================================
    # 7. FINAL GEOMETRY CHECK AND CONTIGUOUS RENUMBERING
    # ========================================================

    check_boundary_topology()

    # --------------------------------------------------------
    # 7.1. Verify the protected vertices
    # --------------------------------------------------------

    for i in protected:

        if (
            i not in points
            or not np.array_equal(points[i], source_points[i])
        ):
            raise RuntimeError(
                'A protected vertex was modified.'
            )

    # --------------------------------------------------------
    # 7.2. Verify the protected boundary edges
    # --------------------------------------------------------

    for e, record in original_boundaries.items():

        if (
            record[2] != 'outlet'
            and boundaries.get(e) != record
        ):
            raise RuntimeError(
                'A protected boundary edge was modified.'
            )

    # --------------------------------------------------------
    # 7.3. Verify physical domain area conservation
    # --------------------------------------------------------

    original_area = sum(
        double_area(*source_points[t])
        for t in source_triangles
    )

    final_area = sum(
        double_area(*coordinates(t, {}))
        for t in triangles.values()
    )

    if not np.isclose(
        original_area,
        final_area,
        rtol=1e-9,
        atol=0
    ):
        raise RuntimeError(
            'Physical domain area changed.'
        )

    # --------------------------------------------------------
    # 7.4. Construct contiguous vertex numbering
    # --------------------------------------------------------

    # Edge collapses introduce gaps in the working vertex
    # numbering.
    #
    # SU2 requires a contiguous zero-based numbering.

    used = sorted({
        i
        for t in triangles.values()
        for i in t
    })

    if set(used) != set(points):
        raise RuntimeError(
            'Unused working vertex detected.'
        )

    ids = {
        i: j
        for j, i in enumerate(used)
    }

    # Renumber coordinates.

    new_points = np.asarray([
        points[i]
        for i in used
    ])

    # Renumber triangular connectivity.

    new_triangles = np.asarray(
        [
            [ids[i] for i in t]
            for t in triangles.values()
        ],
        dtype=np.int64
    )

    # Renumber boundary connectivity.

    new_edges = np.asarray(
        [
            [ids[a], ids[b]]
            for a, b, _ in boundaries.values()
        ],
        dtype=np.int64
    )

    new_tags = np.asarray([
        tag
        for _, _, tag in boundaries.values()
    ])

    # --------------------------------------------------------
    # 7.5. Verify final triangle orientation
    # --------------------------------------------------------

    if any(
        double_area(*new_points[t]) <= 0
        for t in new_triangles
    ):
        raise RuntimeError(
            'Invalid final triangle orientation.'
        )

    # --------------------------------------------------------
    # 7.6. Compute physical triangle quality
    # --------------------------------------------------------

    # Report Euclidean quality separately:
    #
    # Metric-aligned triangles can be highly stretched in
    # physical coordinates.

    x = new_points[new_triangles]

    u = x[:, 1] - x[:, 0]

    v = x[:, 2] - x[:, 1]

    w = x[:, 0] - x[:, 2]

    area2 = (
        u[:, 0] * (-w[:, 1])
        - u[:, 1] * (-w[:, 0])
    )

    physical_quality = (
        2
        * np.sqrt(3)
        * area2
        / np.sum(u*u + v*v + w*w, axis=1)
    )

    # Identify frozen edges that cannot satisfy the
    # prescribed metric length interval.

    frozen_violations = sum(
        not collapse_threshold <= L(e) <= split_threshold
        for e, record in boundaries.items()
        if record[2] != 'outlet'
    )

    print(
        f'Physical quality: '
        f'min={physical_quality.min():.6g}, '
        f'median={np.median(physical_quality):.6g}',
        flush=True
    )

    print(
        f'Frozen boundary edges outside metric thresholds: '
        f'{frozen_violations}',
        flush=True
    )

    print(
        f'Remeshing status: {status}',
        flush=True
    )

    if diagnostics is not None:
        diagnostics.update(totals)
        diagnostics.update(
            metric_quality_min=qmin,
            metric_quality_median=qmed,
            physical_quality_min=float(physical_quality.min()),
            physical_quality_median=float(np.median(physical_quality)),
            metric_edge_rms=error,
            short_edge_count=short,
            long_edge_count=long,
            frozen_boundary_violations=frozen_violations,
        )

    return (
        new_points,
        new_triangles,
        new_edges,
        new_tags,
        status
    )


# ============================================================
# 8. SU2 EXPORT
# ============================================================

def write_su2(
    output_path,
    points,
    triangles,
    boundary_edges,
    boundary_tags
):
    """
    Write native SU2 ASCII connectivity and marker groups
    directly.

    No Gmsh recompilation or external mesh generator is
    required.
    """

    output_path = Path(output_path)

    if output_path.suffix.lower() != '.su2':
        raise ValueError(
            'Output must have a .su2 extension.'
        )

    with output_path.open(
        'w',
        encoding='utf-8'
    ) as f:

        # ----------------------------------------------------
        # Spatial dimension and volume element connectivity
        # ----------------------------------------------------

        f.write(
            f'NDIME= 2\n'
            f'NELEM= {len(triangles)}\n'
        )

        for T, (i, j, l) in enumerate(triangles):

            f.write(
                f'5 {i} {j} {l} {T}\n'
            )

        # ----------------------------------------------------
        # Vertex coordinates
        # ----------------------------------------------------

        f.write(
            f'NPOIN= {len(points)}\n'
        )

        for i, (x, y) in enumerate(points):

            f.write(
                f'{x:.17g} {y:.17g} {i}\n'
            )

        # ----------------------------------------------------
        # Boundary marker definitions
        # ----------------------------------------------------

        tags = list(
            dict.fromkeys(boundary_tags)
        )

        f.write(
            f'NMARK= {len(tags)}\n'
        )

        for tag in tags:

            edges = boundary_edges[
                boundary_tags == tag
            ]

            f.write(
                f'MARKER_TAG= {tag}\n'
                f'MARKER_ELEMS= {len(edges)}\n'
            )

            for a, b in edges:

                f.write(
                    f'3 {a} {b}\n'
                )


# ============================================================
# 9. MAIN PROGRAM
# ============================================================

def mesh_adaptation(
    input_path,
    target_nodes,
    epsilon,
    p,
    h_min,
    h_max,
    max_aspect_ratio,
    collapse_threshold,
    split_threshold,
    max_iterations,
    output_path='adapted_mesh.su2',
    *,
    diagnostics=None
):
    """
    Construct the optimal bounded metric and adapt the mesh.

    Parameters
    ----------
    input_path : str or Path
        NPZ file containing the original mesh and scalar field.

    target_nodes : int
        Indicative number of vertices used to calculate the
        target metric complexity.

    epsilon : float
        Strictly positive Hessian curvature floor.

    p : int
        Exponent of the Lebesgue Lp norm.

    h_min : float
        Minimum prescribed physical metric size.

    h_max : float
        Maximum prescribed physical metric size.

    max_aspect_ratio : float
        Maximum ratio between the principal metric sizes.

    collapse_threshold : float
        Metric length below which an edge can be contracted.

    split_threshold : float
        Metric length above which an edge can be subdivided.

    max_iterations : int
        Maximum number of remeshing passes.

    output_path : str or Path
        Destination SU2 mesh file.

    Returns
    -------
    status : str
        converged, stationary or iteration_limit.

    Notes
    -----
    The optional diagnostics dictionary receives already computed statistics.
    The return value and mathematical procedure are unchanged.
    """

    # --------------------------------------------------------
    # 9.1. Validate adaptation parameters
    # --------------------------------------------------------

    if (
        not isinstance(target_nodes, (int, np.integer))
        or target_nodes <= 0
        or not isinstance(p, (int, np.integer))
        or p <= 0
        or not isinstance(max_iterations, (int, np.integer))
        or max_iterations <= 0
    ):
        raise ValueError(
            'target_nodes, p and max_iterations '
            'must be positive integers.'
        )

    if not np.isfinite([
        epsilon,
        h_min,
        h_max,
        max_aspect_ratio,
        collapse_threshold,
        split_threshold
    ]).all():

        raise ValueError(
            'Non-finite adaptation parameter.'
        )

    if not (
        epsilon > 0
        and 0 < h_min <= h_max
        and max_aspect_ratio >= 1
        and 0 < collapse_threshold < 1 < split_threshold
    ):
        raise ValueError(
            'Invalid metric or remeshing parameters.'
        )

    if Path(input_path).resolve() == Path(output_path).resolve():
        raise ValueError(
            'Output would overwrite input.'
        )

    # --------------------------------------------------------
    # 9.2. Load the source mesh and scalar field
    # --------------------------------------------------------

    (
        points,
        triangles,
        scalar,
        edges,
        tags,
        protected,
        field
    ) = load_input(input_path)

    # --------------------------------------------------------
    # 9.3. Compute the source triangle geometry
    # --------------------------------------------------------

    B, areas = triangle_geometry(
        points,
        triangles
    )

    incident = incident_triangles(
        triangles,
        len(points)
    )

    # Compute the sum of the physical areas of all
    # triangles incident to each vertex.

    area_sums = np.asarray([
        areas[ids].sum()
        for ids in incident
    ])

    if np.any(area_sums <= 0):
        raise ValueError(
            'Isolated source vertex.'
        )

    # Lumped nodal areas:
    #
    # A_i = (1/3) * sum_{T containing i} |T|

    A = area_sums / 3.0

    # --------------------------------------------------------
    # 9.4. Double P1 Hessian reconstruction
    # --------------------------------------------------------

    H = recover_hessian(
        triangles,
        scalar,
        B,
        incident,
        areas,
        area_sums
    )

    # --------------------------------------------------------
    # 9.5. Regularised curvature tensors
    # --------------------------------------------------------

    K, det_K = regularised_curvature(
        H,
        epsilon
    )

    # Release intermediate arrays that are no longer
    # required after the curvature reconstruction.

    del H, scalar, B, incident, area_sums

    # --------------------------------------------------------
    # 9.6. Analytical optimal metric
    # --------------------------------------------------------

    M_star, C_star, s_star = optimal_metric(
        A,
        K,
        det_K,
        target_nodes,
        p
    )

    del K, det_K

    C_optimal = complexity(
        A,
        M_star
    )

    # --------------------------------------------------------
    # 9.7. Metric bounding and complexity restoration
    # --------------------------------------------------------

    M_b, beta = bounded_metric(
        M_star,
        A,
        C_star,
        h_min,
        h_max,
        max_aspect_ratio
    )

    del M_star

    C_bounded = complexity(
        A,
        M_b
    )

    if diagnostics is not None:
        diagnostics.update(
            bounded_metric_complexity=C_bounded,
            complexity_relative_error=abs(C_bounded-C_star)/C_star,
        )

    # Verify that both the analytical optimal metric
    # and the bounded metric respect the target complexity.

    if (
        abs(C_optimal-C_star)/C_star > 1e-8
        or abs(C_bounded-C_star)/C_star > 1e-8
    ):
        raise RuntimeError(
            'Metric complexity verification failed.'
        )

    print(
        f'Field: {field}\n'
        f'Target complexity: {C_star:.12g}\n'
        f'Optimal complexity: {C_optimal:.12g}\n'
        f'Bounded complexity: {C_bounded:.12g}\n'
        f's_star: {s_star:.6g}, '
        f'beta: {beta:.6g}',
        flush=True
    )

    # --------------------------------------------------------
    # 9.8. Construct the immutable source metric field
    # --------------------------------------------------------

    (
        metric_at,
        metric_length,
        metric_midpoint
    ) = build_metric_field(
        points,
        triangles,
        M_b
    )

    # --------------------------------------------------------
    # 9.9. Anisotropic remeshing
    # --------------------------------------------------------

    result = adapt_mesh(
        points,
        triangles,
        edges,
        tags,
        protected,
        metric_at,
        metric_length,
        metric_midpoint,
        collapse_threshold,
        split_threshold,
        max_iterations,
        diagnostics=diagnostics
    )

    (
        new_points,
        new_triangles,
        new_edges,
        new_tags,
        status
    ) = result

    # --------------------------------------------------------
    # 9.10. Export the adapted mesh
    # --------------------------------------------------------

    write_su2(
        output_path,
        new_points,
        new_triangles,
        new_edges,
        new_tags
    )

    print(
        f'Final nodes: {len(new_points)} '
        f'(indicative target: {target_nodes})\n'
        f'Final triangles: {len(new_triangles)}\n'
        f'Output: {output_path}',
        flush=True
    )

    return status


# ============================================================
# 10. COMMAND-LINE INTERFACE
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        'input_path'
    )

    parser.add_argument(
        'target_nodes',
        type=int
    )

    parser.add_argument(
        'epsilon',
        type=float
    )

    parser.add_argument(
        'p',
        type=int
    )

    parser.add_argument(
        'h_min',
        type=float
    )

    parser.add_argument(
        'h_max',
        type=float
    )

    parser.add_argument(
        'max_aspect_ratio',
        type=float
    )

    parser.add_argument(
        'collapse_threshold',
        type=float
    )

    parser.add_argument(
        'split_threshold',
        type=float
    )

    parser.add_argument(
        'max_iterations',
        type=int
    )

    parser.add_argument(
        'output_path',
        nargs='?',
        default='adapted_mesh.su2'
    )

    args = parser.parse_args()

    try:

        mesh_adaptation(
            **vars(args)
        )

    except (
        ValueError,
        RuntimeError,
        OSError,
        KeyError,
        np.linalg.LinAlgError
    ) as error:

        parser.exit(
            1,
            f'Error: {error}\n'
        )


if __name__ == '__main__':
    main()
