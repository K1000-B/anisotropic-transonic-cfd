"""Immutable source mesh and connectivity for a conforming planar triangulation."""

from collections import defaultdict
from dataclasses import dataclass, field
from types import MappingProxyType
from collections.abc import Mapping

import numpy as np


Edge = tuple[int, int]
ROUND_OFF = 64 * np.finfo(float).eps


def edge_key(a: int, b: int) -> Edge:
    return (a, b) if a < b else (b, a)


def triangle_edges(triangle):
    a, b, c = triangle
    return (edge_key(a, b), edge_key(b, c), edge_key(c, a))


def cross(u, v):
    """Signed 2D cross product with NumPy broadcasting."""
    return u[..., 0] * v[..., 1] - u[..., 1] * v[..., 0]


def positive_triangles(points):
    """Check CCW orientation above a local floating-point area tolerance."""
    u = points[..., 1, :] - points[..., 0, :]
    v = points[..., 2, :] - points[..., 0, :]
    longest_squared = np.maximum.reduce(
        (np.sum(u * u, axis=-1), np.sum(v * v, axis=-1), np.sum((u - v) ** 2, axis=-1))
    )
    signed_area = cross(u, v)
    return np.isfinite(signed_area) & (signed_area > ROUND_OFF * longest_squared)


def readonly(values):
    data = np.ascontiguousarray(values)
    return np.frombuffer(data.tobytes(), dtype=data.dtype).reshape(data.shape)


def integer_array(values, label):
    array = np.asarray(values)
    if array.dtype.kind not in "iu":
        raise ValueError(f"{label} must contain integer indices")
    if array.size and (array.min() < 0 or array.max() > np.iinfo(np.int64).max):
        raise ValueError(f"{label} contains an invalid index")
    return array.astype(np.int64, copy=False)


def reachable(start, neighbors):
    visited = {start}
    pending = [start]
    while pending:
        for item in neighbors[pending.pop()]:
            if item not in visited:
                visited.add(item)
                pending.append(item)
    return visited


@dataclass(frozen=True, eq=False)
class Mesh:
    """Arrays (N,2), (T,3), and mapping of boundary names to (E_marker,2).

    Construct with counterclockwise triangles and exactly all boundary edges.
    Point IDs are optional external identifiers, independent of array indices.
    """

    coordinates: np.ndarray
    triangles: np.ndarray
    markers: Mapping[str, np.ndarray]
    point_ids: np.ndarray | None = None
    triangle_areas: np.ndarray = field(init=False)
    nodal_areas: np.ndarray = field(init=False)
    incidence: tuple = field(init=False)
    neighbors: tuple = field(init=False)
    edges: Mapping = field(init=False)
    boundary: Mapping = field(init=False)
    boundary_components: int = field(init=False)
    components: int = field(init=False)
    euler_characteristic: int = field(init=False)

    def __post_init__(self):
        xy = np.asarray(self.coordinates, dtype=float)
        tri = integer_array(self.triangles, "triangles")
        if xy.ndim != 2 or xy.shape[1] != 2 or not np.isfinite(xy).all():
            raise ValueError("coordinates must be finite with shape (N, 2)")
        n = len(xy)
        if tri.ndim != 2 or tri.shape[1] != 3 or len(tri) == 0 or tri.max() >= n:
            raise ValueError("invalid triangular connectivity")
        if len(np.unique(xy, axis=0)) != n:
            raise ValueError("duplicate vertex coordinates")
        if len(np.unique(np.sort(tri, axis=1), axis=0)) != len(tri):
            raise ValueError("duplicate triangle")
        points = xy[tri]
        if not positive_triangles(points).all():
            raise ValueError("triangles must be nondegenerate and counterclockwise")
        x0, x1, x2 = points.transpose(1, 0, 2)
        areas = cross(x1 - x0, x2 - x0) / 2
        lumped = np.bincount(tri.ravel(), np.repeat(areas / 3, 3), minlength=n)
        if np.any(lumped <= 0) or not np.isfinite(lumped.sum()):
            raise ValueError("unused vertex or invalid total area")

        edges = defaultdict(list)
        incidence = [[] for _ in xy]
        neighbors = [set() for _ in xy]
        for i, (a, b, c) in enumerate(tri.tolist()):
            for vertex in (a, b, c):
                incidence[vertex].append(i)
            for start, end in ((a, b), (b, c), (c, a)):
                edges[edge_key(start, end)].append((i, start < end))
                neighbors[start].add(end)
                neighbors[end].add(start)
        if any(len(owners) > 2 or (len(owners) == 2 and owners[0][1] == owners[1][1])
               for owners in edges.values()):
            raise ValueError("non-manifold or inconsistently oriented edge")

        markers = {}
        boundary = {}
        rim = defaultdict(set)
        for marker, values in self.markers.items():
            if not isinstance(marker, str) or not marker or marker.strip() != marker or "\n" in marker:
                raise ValueError("invalid boundary marker name")
            segments = integer_array(values, "boundary segments")
            if segments.ndim != 2 or segments.shape[1] != 2 or not len(segments):
                raise ValueError("marker connectivity must have shape (E, 2), E > 0")
            for a, b in segments.tolist():
                edge = edge_key(a, b)
                if edge in boundary or len(edges.get(edge, ())) != 1:
                    raise ValueError(f"invalid or repeated boundary edge {edge}")
                boundary[edge] = marker
                rim[a].add(b)
                rim[b].add(a)
            markers[marker] = readonly(segments)
        if len(boundary) != sum(len(owners) == 1 for owners in edges.values()):
            raise ValueError("markers must cover all and only one-sided edges")
        if any(len(adjacent) != 2 for adjacent in rim.values()):
            raise ValueError("boundary curves must be disjoint loops in connectivity")

        for vertex, faces in enumerate(incidence):
            link = defaultdict(set)
            for face in faces:
                a, b = (v for v in tri[face].tolist() if v != vertex)
                link[a].add(b)
                link[b].add(a)
            if len(reachable(next(iter(link)), link)) != len(link):
                raise ValueError(f"disconnected triangle fan at vertex {vertex}")

        self._check_boundary_geometry(xy, boundary)
        self._check_boundary_winding(xy, boundary, edges)
        ids = np.arange(n) if self.point_ids is None else integer_array(self.point_ids, "point IDs")
        if ids.shape != (n,) or len(set(ids.tolist())) != n:
            raise ValueError("point IDs must be one per vertex and unique")
        attributes = dict(
            coordinates=readonly(xy), triangles=readonly(tri), markers=MappingProxyType(markers),
            point_ids=readonly(ids), triangle_areas=readonly(areas), nodal_areas=readonly(lumped),
            incidence=tuple(tuple(faces) for faces in incidence),
            neighbors=tuple(frozenset(adjacent) for adjacent in neighbors),
            edges=MappingProxyType({e: tuple(i for i, _ in owners) for e, owners in edges.items()}),
            boundary=MappingProxyType(boundary), components=self._count_components(range(n), neighbors),
            boundary_components=self._count_components(rim, rim),
            euler_characteristic=n - len(edges) + len(tri),
        )
        for name, value in attributes.items():
            object.__setattr__(self, name, value)

    @staticmethod
    def _count_components(vertices, neighbors):
        unseen = set(vertices)
        count = 0
        while unseen:
            unseen.difference_update(reachable(next(iter(unseen)), neighbors))
            count += 1
        return count

    @staticmethod
    def _check_boundary_geometry(xy, boundary):
        """Sweep boundary boxes and reject crossings and collinear overlaps."""
        segments = np.array(sorted(boundary), dtype=int).reshape(-1, 2)
        a, b = xy[segments].transpose(1, 0, 2)
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        order = np.argsort(lo[:, 0], kind="stable")
        active = []
        for i in order:
            active = [j for j in active if hi[j, 0] >= lo[i, 0]]
            for j in active:
                if hi[j, 1] < lo[i, 1] or hi[i, 1] < lo[j, 1]:
                    continue
                shared = set(segments[i]).intersection(segments[j])
                if shared:
                    vertex, = shared
                    u = xy[next(v for v in segments[i] if v != vertex)] - xy[vertex]
                    v = xy[next(v for v in segments[j] if v != vertex)] - xy[vertex]
                    if cross(u, v) == 0 and u @ v > 0:
                        raise ValueError("overlapping incident boundary segments")
                    continue
                p = cross(b[i] - a[i], a[j] - a[i])
                q = cross(b[i] - a[i], b[j] - a[i])
                r = cross(b[j] - a[j], a[i] - a[j])
                s = cross(b[j] - a[j], b[i] - a[j])
                if min(p, q) <= 0 <= max(p, q) and min(r, s) <= 0 <= max(r, s):
                    raise ValueError("intersecting boundary segments")
            active.append(i)

    @staticmethod
    def _check_boundary_winding(xy, boundary, edges):
        """Require winding numbers zero or one on a mesh with simple boundary loops."""
        successor = {}
        for a, b in boundary:
            if edges[a, b][0][1]:
                successor[a] = b
            else:
                successor[b] = a
        unseen = set(successor)
        loops = []
        while unseen:
            start = min(unseen)
            vertices = [start]
            vertex = successor[start]
            while vertex != start:
                vertices.append(vertex)
                vertex = successor[vertex]
            unseen.difference_update(vertices)
            a = xy[vertices]
            b = np.roll(a, -1, axis=0)
            signed_area = cross(a - a[0], b - a[0]).sum()
            if signed_area == 0 or not np.isfinite(signed_area):
                raise ValueError("invalid boundary loop area")
            loops.append((a, b, 1 if signed_area > 0 else -1))
        for i, (loop, _, sign) in enumerate(loops):
            x, y = loop[0]
            winding = 0
            for j, (a, b, orientation) in enumerate(loops):
                if i == j:
                    continue
                crosses_y = (a[:, 1] > y) != (b[:, 1] > y)
                p, q = a[crosses_y], b[crosses_y]
                to_right = p[:, 0] + (y - p[:, 1]) * (q[:, 0] - p[:, 0]) / (q[:, 1] - p[:, 1]) > x
                winding += orientation * (np.count_nonzero(to_right) % 2)
            if winding != (0 if sign > 0 else 1):
                raise ValueError("overlapping domains or invalid boundary winding")
