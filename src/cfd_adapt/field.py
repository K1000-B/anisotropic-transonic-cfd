"""P1 interpolation and piecewise analytic line integration on a fixed source."""

from dataclasses import dataclass

import numpy as np

from .mesh import Mesh, readonly


class OutsideSourceError(ValueError):
    """A point or segment lies outside the polygonal source domain."""


@dataclass(frozen=True)
class _Box:
    lower: tuple
    upper: tuple
    triangles: np.ndarray | None
    children: tuple = ()


class P1MetricField:
    """Evaluate (N,2,2) nodal SPD tensors on an immutable Mesh.

    A bounding volume hierarchy indexes triangle boxes. The dimensionless
    tolerance admits round-off at triangle boundaries. Accepted barycentric
    weights are clipped back to the closed triangle before evaluation.
    """

    def __init__(self, mesh: Mesh, nodal_metric, *, tolerance=1e-10):
        matrices = np.asarray(nodal_metric, dtype=float)
        if (matrices.shape != (len(mesh.coordinates), 2, 2)
                or not np.isfinite(matrices).all()
                or not np.array_equal(matrices, matrices.swapaxes(1, 2))
                or np.any(np.linalg.eigvalsh(matrices) <= 0)):
            raise ValueError("expected one finite symmetric positive definite tensor per source node")
        if not 0 < tolerance < 1e-4:
            raise ValueError("location tolerance must be positive and less than 1e-4")
        self._mesh = mesh
        self._nodal = readonly(matrices)
        self._tolerance = tolerance
        points = mesh.coordinates[mesh.triangles]
        self._origins = points[:, 0]
        self._inverse = np.linalg.inv((points[:, 1:] - points[:, :1]).transpose(0, 2, 1))
        lo, hi = points.min(axis=1), points.max(axis=1)
        padding = 2 * tolerance * (hi - lo)
        self._tree = self._build(np.arange(len(points)), lo - padding, hi + padding)

    @property
    def mesh(self):
        return self._mesh

    @property
    def nodal(self):
        return self._nodal

    @property
    def tolerance(self):
        return self._tolerance

    def _build(self, indices, lower, upper):
        lo, hi = lower[indices].min(axis=0), upper[indices].max(axis=0)
        if len(indices) <= 8:
            return _Box(tuple(lo), tuple(hi), indices)
        axis = int(np.argmax(hi - lo))
        order = np.argsort(lower[indices, axis] + upper[indices, axis], kind="stable")
        left, right = np.split(indices[order], [len(indices) // 2])
        return _Box(tuple(lo), tuple(hi), None,
                    (self._build(left, lower, upper), self._build(right, lower, upper)))

    def _candidates(self, a, b=None):
        x, y = map(float, a)
        dx, dy = (0.0, 0.0) if b is None else map(float, b - a)
        inv_x = 1.0 / dx if dx else 0.0
        inv_y = 1.0 / dy if dy else 0.0
        found = []
        pending = [self._tree]
        while pending:
            box = pending.pop()
            x0, y0 = box.lower
            x1, y1 = box.upper
            near, far = 0.0, 1.0
            if dx:
                x_enter, x_exit = (x0, x1) if dx > 0 else (x1, x0)
                near, far = max(near, (x_enter - x) * inv_x), min(far, (x_exit - x) * inv_x)
            elif x < x0 or x > x1:
                continue
            if dy:
                y_enter, y_exit = (y0, y1) if dy > 0 else (y1, y0)
                near, far = max(near, (y_enter - y) * inv_y), min(far, (y_exit - y) * inv_y)
            elif y < y0 or y > y1:
                continue
            if near <= far:
                if box.triangles is None:
                    pending.extend(box.children)
                else:
                    found.append(box.triangles)
        return np.concatenate(found) if found else np.empty(0, dtype=int)

    @staticmethod
    def _barycentric(tail):
        return np.column_stack((1 - tail.sum(axis=-1), tail))

    @staticmethod
    def _point(values):
        point = np.asarray(values, dtype=float)
        if point.shape != (2,) or not np.isfinite(point).all():
            raise ValueError("expected a finite planar coordinate")
        return point

    @staticmethod
    def _closed_weights(values):
        weights = np.maximum(values, 0)
        return weights / weights.sum(axis=-1, keepdims=True)

    def locate(self, point):
        """Return a source triangle index and three closed barycentric weights."""
        point = self._point(point)
        candidates = self._candidates(point)
        tail = np.einsum("tij,tj->ti", self._inverse[candidates], point - self._origins[candidates])
        weights = self._barycentric(tail)
        lower = weights.min(axis=1)
        if len(lower):
            best = int(np.argmax(lower))
            if lower[best] >= -self.tolerance:
                return int(candidates[best]), self._closed_weights(weights[best])
        raise OutsideSourceError(f"point {point.tolist()} lies outside the source domain")

    def at(self, point):
        """eq:interpolated_metric: evaluate the barycentric P1 tensor at (2,)."""
        face, weights = self.locate(point)
        return np.einsum("i,ijk->jk", weights, self.nodal[self.mesh.triangles[face]])

    def _segments(self, a, b):
        """Clip the line in barycentric space and return nonoverlapping intervals."""
        candidates = self._candidates(a, b)
        inverse = self._inverse[candidates]
        origins = self._origins[candidates]
        w_a = self._barycentric(np.einsum("tij,tj->ti", inverse, a - origins))
        w_b = self._barycentric(np.einsum("tij,tj->ti", inverse, b - origins))
        close_to_side = np.maximum(abs(w_a), abs(w_b)) <= self.tolerance
        delta = np.where(close_to_side, 0, w_b - w_a)
        values = np.where(close_to_side, 0, w_a)
        roots = np.divide(-values, delta, out=np.zeros_like(values), where=delta != 0)
        lower = np.maximum(0, np.max(np.where(delta > 0, roots, -np.inf), axis=1))
        upper = np.minimum(1, np.min(np.where(delta < 0, roots, np.inf), axis=1))
        valid = (upper > lower) & ~np.any((delta == 0) & (values < 0), axis=1)
        order = np.flatnonzero(valid)
        order = order[np.lexsort((-upper[order], lower[order]))]
        end = 0.0
        segments = []
        for index in order:
            if upper[index] <= end:
                continue
            if lower[index] - end > self.tolerance:
                break
            stop = min(1.0, float(upper[index]))
            if 1.0 - stop <= self.tolerance:
                stop = 1.0
            values = np.stack((w_a[index], w_b[index]))
            weights = self._closed_weights(values[0] + np.array([[end], [stop]]) * (values[1] - values[0]))
            segments.append((end, stop, candidates[index], weights))
            end = stop
            if end == 1.0:
                return segments
        raise OutsideSourceError(f"segment {a.tolist()} to {b.tolist()} leaves the source domain")

    @staticmethod
    def _linear_sqrt_mean(first, last):
        # Let a and b be endpoint roots. Cancelling (b-a) gives a stable
        # integral even as both endpoint values approach one another.
        lo, hi = np.minimum(first, last), np.maximum(first, last)
        ratio = lo / hi
        return (2 / 3) * hi * (1 + ratio * ratio / (1 + ratio))

    @staticmethod
    def _directional_roots(matrices, direction):
        scale = np.max(np.abs(direction))
        unit = direction / scale
        squared = np.einsum("...ij,i,j->...", matrices, unit, unit)
        if np.any(squared <= 0) or not np.isfinite(squared).all():
            raise ArithmeticError("invalid directional metric in edge integration")
        return np.sqrt(squared) * scale

    def source_edge_lengths(self, edges):
        """Vectorised analytic lengths for source edges or subsegments of one face."""
        edges = np.asarray(edges, dtype=int).reshape(-1, 2)
        directions = self.mesh.coordinates[edges[:, 1]] - self.mesh.coordinates[edges[:, 0]]
        scale = np.max(np.abs(directions), axis=1)
        unit = directions / scale[:, None]
        squared = np.einsum("ekij,ei,ej->ek", self.nodal[edges], unit, unit)
        if np.any(squared <= 0) or not np.isfinite(squared).all():
            raise ArithmeticError("invalid directional metric on source edge")
        roots = np.sqrt(squared) * scale[:, None]
        lengths = self._linear_sqrt_mean(roots[:, 0], roots[:, 1])
        if not np.isfinite(lengths).all():
            raise ArithmeticError("source edge lengths exceed floating-point range")
        return lengths

    def edge_measure(self, a, b, *, midpoint=False):
        """eq:edge_metric_length: integrate every P1 segment analytically.

        For midpoint=True, also return the physical point at half metric length.
        Reversal uses the same canonical orientation and gives identical results.
        """
        a, b = self._point(a), self._point(b)
        if tuple(b) < tuple(a):
            a, b = b, a
        if np.array_equal(a, b):
            self.locate(a)
            return (0.0, a.copy()) if midpoint else 0.0
        segments = self._segments(a, b)
        weights = np.stack([segment[3] for segment in segments])
        faces = np.array([segment[2] for segment in segments])
        endpoints = np.array([segment[:2] for segment in segments])
        metric = np.einsum("skv,svij->skij", weights, self.nodal[self.mesh.triangles[faces]])
        roots = self._directional_roots(metric, b - a)
        widths = endpoints[:, 1] - endpoints[:, 0]
        lengths = widths * self._linear_sqrt_mean(roots[:, 0], roots[:, 1])
        total = float(lengths.sum())
        if not np.isfinite(total) or total <= 0:
            raise ArithmeticError("edge length exceeds floating-point range")
        if not midpoint:
            return total
        ends = np.cumsum(lengths)
        face = int(np.searchsorted(ends, 0.5 * total))
        fraction = (0.5 * total - (ends[face - 1] if face else 0.0)) / lengths[face]
        # On this segment invert the integral using dimensionless endpoint roots.
        first, last = roots[face] / roots[face].max()
        whole = self._linear_sqrt_mean(first, last)
        low, high = 0.0, 1.0
        for _ in range(40):
            t = 0.5 * (low + high)
            root = np.hypot(np.sqrt(1 - t) * first, np.sqrt(t) * last)
            if t * self._linear_sqrt_mean(first, root) < fraction * whole:
                low = t
            else:
                high = t
        position = endpoints[face, 0] + 0.5 * (low + high) * widths[face]
        return total, a + position * (b - a)
