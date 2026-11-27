"""Local split, collapse, swap and relocation on a fixed P1 metric."""

from collections import Counter
from dataclasses import dataclass, field
from numbers import Integral

import numpy as np

from .field import OutsideSourceError, P1MetricField
from .mesh import Mesh, cross, edge_key, positive_triangles, triangle_edges


@dataclass(frozen=True)
class RemeshParameters:
    collapse_below: float = 0.62
    split_above: float = 1.45
    min_quality: float = 0.08
    convergence_quality: float = 0.20
    quality_gain: float = 1e-5
    relocation_steps: int = 4
    max_passes: int = 8
    max_operations: int = 100_000
    location_tolerance: float = 1e-10

    def __post_init__(self):
        numbers = (self.collapse_below, self.split_above, self.min_quality,
                   self.convergence_quality, self.quality_gain, self.location_tolerance)
        if not np.isfinite(numbers).all() or not (
            0 < self.collapse_below < 1 < self.split_above
            and 0 < self.min_quality <= self.convergence_quality <= 1
            and 0 < self.quality_gain < 1 and 0 < self.location_tolerance < 1e-4
        ):
            raise ValueError("invalid remeshing thresholds or tolerances")
        for value in (self.relocation_steps, self.max_passes, self.max_operations):
            if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
                raise ValueError("remeshing step, pass and operation limits must be positive integers")


@dataclass(frozen=True)
class GeometryConstraints:
    """Input vertex indices and edges that must retain their geometry.

    By default every boundary edge is held fixed. An explicit allow-list
    permits splitting straight edges. None permits splitting the full input
    polygon, including chords of curved CAD boundaries if present.
    """

    fixed_vertices: frozenset[int] = frozenset()
    protected_edges: frozenset[tuple[int, int]] = frozenset()
    splittable_boundary_edges: frozenset[tuple[int, int]] | None = frozenset()


@dataclass
class RemeshDiagnostics:
    termination: str = ""
    passes: list[dict] = field(default_factory=list)
    accepted: Counter = field(default_factory=Counter)
    rejected: Counter = field(default_factory=Counter)
    initial_edge_lengths: np.ndarray | None = None
    final_edge_lengths: np.ndarray | None = None
    initial_edges: np.ndarray | None = None
    final_edges: np.ndarray | None = None
    initial_physical_quality: np.ndarray | None = None
    final_physical_quality: np.ndarray | None = None
    initial_metric_quality: np.ndarray | None = None
    final_metric_quality: np.ndarray | None = None
    initial_outside: tuple[int, int] = (0, 0)
    final_outside: tuple[int, int] = (0, 0)
    initial_edge_penalty: float = 0.0
    final_edge_penalty: float = 0.0
    final_constrained_edges: int = 0


def triangle_quality(points, metric=None):
    """MA-T01: mean-ratio quality for physical or constant-metric (T,3,2) triangles.

    Metric quality uses the mean of the three sampled vertex tensors. This
    auxiliary frozen tensor is independent of the exact P1 edge integral.
    """
    a, b, c = np.moveaxis(np.asarray(points), -2, 0)
    u, v = b - a, c - a
    sides = np.stack((u, v, v - u), axis=-2)
    area2 = np.abs(cross(u, v))
    if metric is None:
        squared = np.sum(sides * sides, axis=(-2, -1))
    else:
        scale = np.max(np.abs(metric), axis=(-2, -1))
        normalized = metric / scale[..., None, None]
        squared = np.einsum("...ki,...ij,...kj->...", sides, normalized, sides)
        area2 = area2 * np.sqrt(np.linalg.det(normalized))
    return 2 * np.sqrt(3) * area2 / squared


class _Triangulation:
    """Maintain incident triangles, edge owners and caches through local edits."""

    def __init__(self, field, constraints):
        mesh = field.mesh
        self.field = field
        self.xy = dict(enumerate(mesh.coordinates))
        self.metric = dict(enumerate(field.nodal))
        self.triangles = {}
        self.incidence = {v: set() for v in self.xy}
        self.neighbors = {v: set() for v in self.xy}
        self.edges = {}
        self.boundary = dict(mesh.boundary)
        self.boundary_vertices = set(v for edge in self.boundary for v in edge)
        self.protected = {edge_key(*edge) for edge in constraints.protected_edges}
        fixed = set(constraints.fixed_vertices)
        if not self.protected <= mesh.edges.keys() or not fixed <= self.xy.keys():
            raise ValueError("protected edges and fixed vertices must exist in the source mesh")
        self.fixed = fixed | {v for edge in self.protected for v in edge}
        if constraints.splittable_boundary_edges is None:
            self.splittable = set(self.boundary)
        else:
            self.splittable = {edge_key(*edge) for edge in constraints.splittable_boundary_edges}
            if not self.splittable <= self.boundary.keys():
                raise ValueError("only existing boundary edges can be marked splittable")
        self.splittable.difference_update(self.protected)
        self.next_face = 0
        self.next_vertex = len(self.xy)
        self.add_faces(map(tuple, mesh.triangles.tolist()))
        self.lengths = dict(zip(self.edges, field.source_edge_lengths(list(self.edges))))
        self.born = set()

    def add_faces(self, faces):
        for tri in faces:
            face = self.next_face
            self.next_face += 1
            self.triangles[face] = tri
            for vertex in tri:
                self.incidence[vertex].add(face)
            for a, b in triangle_edges(tri):
                if (a, b) not in self.edges:
                    self.edges[a, b] = set()
                    self.neighbors[a].add(b)
                    self.neighbors[b].add(a)
                self.edges[a, b].add(face)

    def replace(self, old, new):
        old = sorted(old)
        touched = set()
        for face in old:
            tri = self.triangles.pop(face)
            for vertex in tri:
                self.incidence[vertex].remove(face)
            for edge in triangle_edges(tri):
                self.edges[edge].remove(face)
                touched.add(edge)
        # Keep caches and adjacencies for edges that survive the local edit.
        self.add_faces(new)
        for a, b in touched:
            if not self.edges[a, b]:
                del self.edges[a, b]
                self.lengths.pop((a, b), None)
                self.neighbors[a].remove(b)
                self.neighbors[b].remove(a)

    def length(self, edge):
        if edge not in self.lengths:
            a, b = edge
            self.lengths[edge] = self.field.edge_measure(self.xy[a], self.xy[b])
        return self.lengths[edge]

    def geometry(self, triangles, vertex=None, position=None, metric=None):
        points = np.array([[position if v == vertex else self.xy[v] for v in tri] for tri in triangles])
        tensors = np.array([[metric if v == vertex else self.metric[v] for v in tri] for tri in triangles])
        return points, tensors.mean(axis=1)

    def quality(self, triangles):
        return triangle_quality(*self.geometry(triangles))

    def snapshot(self):
        keep = sorted(self.xy)
        index = {v: i for i, v in enumerate(keep)}
        triangles = [[index[v] for v in tri] for tri in self.triangles.values()]
        markers = {name: [] for name in self.field.mesh.markers}
        for (a, b), name in sorted(self.boundary.items()):
            markers[name].append((index[a], index[b]))
        mesh = Mesh([self.xy[v] for v in keep], triangles, markers)
        origin = np.array([v if v < len(self.field.mesh.coordinates) else -1 for v in keep])
        return mesh, origin

    def split(self, edge):
        if edge in self.protected or (edge in self.boundary and edge not in self.splittable):
            return "constrained"
        a, b = edge
        _, position = self.field.edge_measure(self.xy[a], self.xy[b], midpoint=True)
        return self._split_at(edge, position)

    def _split_at(self, edge, position):
        a, b = edge
        old = self.edges[edge].copy()
        vertex = self.next_vertex
        new = []
        for face in sorted(old):
            tri = self.triangles[face]
            start = next(i for i, v in enumerate(tri) if v in edge and tri[(i + 1) % 3] in edge)
            u, v, w = (tri[(start + i) % 3] for i in range(3))
            new.extend(((u, vertex, w), (vertex, v, w)))
        point_metric = self.field.at(position)
        points, _ = self.geometry(new, vertex, position, point_metric)
        if not positive_triangles(points).all():
            return "orientation"
        self.next_vertex += 1
        self.xy[vertex] = position
        self.metric[vertex] = point_metric
        self.incidence[vertex] = set()
        self.neighbors[vertex] = set()
        self.born.add(vertex)
        self.replace(old, new)
        if edge in self.boundary:
            name = self.boundary.pop(edge)
            self.splittable.remove(edge)
            for child in (edge_key(a, vertex), edge_key(b, vertex)):
                self.splittable.add(child)
                self.boundary[child] = name
            self.boundary_vertices.add(vertex)
        return None

    def collapse(self, edge, config):
        a, b = edge
        if edge in self.protected or edge in self.boundary:
            return "constrained"
        if a in self.born or b in self.born:
            return "recent_split"
        owners = self.edges[edge]
        opposite = {v for face in owners for v in self.triangles[face] if v not in edge}
        if self.neighbors[a] & self.neighbors[b] != opposite:
            return "link_condition"
        # The two full vertex links must intersect only in the edge link.
        # In particular, a common opposite edge would create duplicate faces.
        if len(opposite) == 2:
            u, v = sorted(opposite)
            common = [self.triangles[face] for face in self.edges.get((u, v), ())]
            if any(a in tri for tri in common) and any(b in tri for tri in common):
                return "link_condition"

        candidates = []
        rejection = "constrained"
        for keep, remove in ((a, b), (b, a)):
            if remove in self.boundary_vertices or remove in self.fixed:
                continue
            changed = self.incidence[remove]
            surviving = [self.triangles[face] for face in sorted(changed - owners)]
            if not surviving:
                rejection = "link_condition"
                continue
            replacement = [tuple(keep if v == remove else v for v in tri) for tri in surviving]
            points, metrics = self.geometry(replacement)
            if not positive_triangles(points).all():
                rejection = "orientation"
                continue
            quality = triangle_quality(points, metrics)
            minimum = float(quality.min())
            patch = [self.triangles[face] for face in sorted(changed | self.incidence[keep])]
            if minimum < min(config.min_quality, self.quality(patch).min()) * (1 - config.quality_gain):
                rejection = "quality"
                continue
            candidates.append((minimum, float(quality.mean()), keep, remove, replacement))

        for _, _, keep, remove, replacement in sorted(candidates, reverse=True):
            new_lengths = {}
            try:
                for vertex in sorted(self.neighbors[remove] - self.neighbors[keep] - {keep}):
                    new_edge = edge_key(keep, vertex)
                    length = self.field.edge_measure(self.xy[keep], self.xy[vertex])
                    before = self.length(edge_key(remove, vertex))
                    if length > max(config.split_above, before) * (1 + config.quality_gain):
                        rejection = "long_replacement"
                        break
                    new_lengths[new_edge] = length
                else:
                    self.replace(self.incidence[remove].copy(), replacement)
                    self.lengths.update(new_lengths)
                    del self.xy[remove], self.metric[remove], self.incidence[remove], self.neighbors[remove]
                    return None
            except OutsideSourceError:
                rejection = "outside_source"
        return rejection

    def swap(self, edge, config):
        if edge in self.protected or edge in self.boundary:
            return "constrained"
        first, second = sorted(self.edges[edge])
        left, right = self.triangles[first], self.triangles[second]
        a, b = edge
        if (left.index(a) + 1) % 3 != left.index(b):
            a, b = b, a
        c = next(v for v in left if v not in edge)
        d = next(v for v in right if v not in edge)
        diagonal = edge_key(c, d)
        if diagonal in self.edges:
            return "link_condition"
        replacement = [(c, a, d), (d, b, c)]
        points, tensors = self.geometry(replacement)
        if not positive_triangles(points).all():
            return "orientation"
        old_quality = self.quality([left, right])
        if triangle_quality(points, tensors).min() <= old_quality.min() + config.quality_gain:
            return "quality"
        length = self.field.edge_measure(self.xy[c], self.xy[d])
        old_length = self.length(edge)
        if length > max(config.split_above, old_length) or length < min(config.collapse_below, old_length):
            return "diagonal_length"
        self.replace([first, second], replacement)
        self.lengths[diagonal] = length
        return None

    def relocate(self, vertex, config):
        if vertex in self.boundary_vertices or vertex in self.fixed:
            return "constrained"
        neighbors = sorted(self.neighbors[vertex])
        triangles = [self.triangles[face] for face in sorted(self.incidence[vertex])]
        points, metrics = self.geometry(triangles)
        old_quality = triangle_quality(points, metrics)
        old_min = old_quality.min()
        position = self.xy[vertex]
        coordinates = np.array([self.xy[v] for v in neighbors])
        tensors = np.array([self.metric[v] for v in neighbors])
        # Minimise the sum of squared frozen one-edge metric distances.
        displacement = np.linalg.solve(tensors.sum(axis=0),
                                       np.einsum("nij,nj->i", tensors, coordinates - position))
        old_lengths = np.array([self.length(edge_key(vertex, v)) for v in neighbors])
        old_penalty = np.mean(np.log(old_lengths) ** 2)
        rejection = "quality"
        for step in range(config.relocation_steps):
            trial = position + 0.5 ** step * displacement
            if np.array_equal(trial, position):
                break
            try:
                metric = self.field.at(trial)
                points, tensors = self.geometry(triangles, vertex, trial, metric)
                if not positive_triangles(points).all():
                    rejection = "orientation"
                    continue
                quality = triangle_quality(points, tensors)
                if quality.min() < min(config.min_quality, old_min) * (1 - config.quality_gain):
                    rejection = "quality"
                    continue
                if quality.min() + quality.mean() < old_min + old_quality.mean() - config.quality_gain:
                    rejection = "quality"
                    continue
                lengths = np.array([self.field.edge_measure(trial, self.xy[v]) for v in neighbors])
            except OutsideSourceError:
                rejection = "outside_source"
                continue
            if (np.any(lengths > np.maximum(config.split_above, old_lengths) * (1 + config.quality_gain))
                    or np.any(lengths < np.minimum(config.collapse_below, old_lengths) * (1 - config.quality_gain))):
                rejection = "edge_length"
                continue
            gain = quality.min() + quality.mean() - old_min - old_quality.mean()
            gain += old_penalty - np.mean(np.log(lengths) ** 2)
            if gain <= config.quality_gain:
                rejection = "quality"
                continue
            self.xy[vertex] = trial
            self.metric[vertex] = metric
            self.lengths.update((edge_key(vertex, v), length) for v, length in zip(neighbors, lengths))
            return None
        return rejection


def _measure(state, config):
    edges = list(state.edges)
    lengths = np.fromiter((state.length(edge) for edge in edges), dtype=float, count=len(edges))
    short = lengths < config.collapse_below
    long = lengths > config.split_above
    points, metrics = state.geometry(list(state.triangles.values()))
    constrained = sum(
        (is_long and (edge in state.protected or (edge in state.boundary and edge not in state.splittable)))
        or (is_short and (edge in state.protected or edge in state.boundary
                          or all(v in state.fixed or v in state.boundary_vertices for v in edge)))
        for edge, is_short, is_long in zip(edges, short, long)
    )
    return (lengths, triangle_quality(points), triangle_quality(points, metrics),
            (int(short.sum()), int(long.sum())), float(np.mean(np.log(lengths) ** 2)), int(constrained))


def _run_pass(state, config, diagnostics, budget):
    accepted = Counter()

    def attempt(name, operation, item):
        try:
            rejection = operation(item)
        except OutsideSourceError:
            rejection = "outside_source"
        if rejection is None:
            accepted[name] += 1
        else:
            diagnostics.rejected[f"{name}.{rejection}"] += 1

    for name, operation, split in (
        ("collapse", lambda edge: state.collapse(edge, config), False),
        ("split", state.split, True),
    ):
        if split:
            state.born.clear()
        priority = -1 if split else 1
        cutoff = config.split_above if split else config.collapse_below
        candidates = sorted((priority * state.length(edge), edge) for edge in state.edges
                            if priority * state.length(edge) < priority * cutoff)
        for _, edge in candidates:
            if accepted.total() >= budget:
                return accepted
            if edge in state.edges and priority * state.length(edge) < priority * cutoff:
                attempt(name, operation, edge)

    for edge in sorted(state.edges):
        if accepted.total() >= budget:
            return accepted
        if len(state.edges.get(edge, ())) == 2:
            attempt("swap", lambda e: state.swap(e, config), edge)

    for vertex in sorted(state.xy.keys() - state.boundary_vertices - state.fixed):
        if accepted.total() >= budget:
            return accepted
        attempt("relocate", lambda v: state.relocate(v, config), vertex)
    return accepted


def remesh(field: P1MetricField, config=RemeshParameters(), constraints=GeometryConstraints()):
    """MA-T01: adapt the source Mesh using an immutable bounded metric field.

    Return (Mesh, diagnostics, source indices). Source indices have length N_out
    and contain -1 for split vertices. A retained index can have moved.
    """
    state = _Triangulation(field, constraints)
    diagnostics = RemeshDiagnostics()
    diagnostics.initial_edges = np.array(list(state.edges))
    (diagnostics.initial_edge_lengths, diagnostics.initial_physical_quality,
     diagnostics.initial_metric_quality, diagnostics.initial_outside,
     diagnostics.initial_edge_penalty, _) = _measure(state, config)
    diagnostics.termination = "pass_limit"
    for sweep in range(config.max_passes):
        remaining = config.max_operations - diagnostics.accepted.total()
        deferred = diagnostics.rejected["collapse.recent_split"]
        accepted = _run_pass(state, config, diagnostics, remaining)
        diagnostics.accepted.update(accepted)
        diagnostics.passes.append({"pass": sweep + 1, "accepted": dict(accepted), "vertices": len(state.xy)})
        if diagnostics.accepted.total() >= config.max_operations:
            diagnostics.termination = "operation_limit"
            break
        if not accepted and diagnostics.rejected["collapse.recent_split"] == deferred:
            diagnostics.termination = "stationary"
            break
    (diagnostics.final_edge_lengths, diagnostics.final_physical_quality,
     diagnostics.final_metric_quality, diagnostics.final_outside,
     diagnostics.final_edge_penalty, diagnostics.final_constrained_edges) = _measure(state, config)
    if (diagnostics.termination == "stationary" and diagnostics.final_outside == (0, 0)
            and diagnostics.final_metric_quality.min() >= config.convergence_quality):
        diagnostics.termination = "criteria_met"

    mesh, origin = state.snapshot()
    index = {v: i for i, v in enumerate(sorted(state.xy))}
    diagnostics.final_edges = np.array([(index[a], index[b]) for a, b in state.edges])
    initial = field.mesh
    if (mesh.components, mesh.boundary_components, mesh.euler_characteristic) != (
        initial.components, initial.boundary_components, initial.euler_characteristic
    ):
        raise ArithmeticError("remeshing changed the domain topology")
    if not np.isclose(mesh.triangle_areas.sum(), initial.triangle_areas.sum(), rtol=2e-13, atol=0):
        raise ArithmeticError("remeshing changed the domain area")
    return mesh, diagnostics, origin
