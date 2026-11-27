"""Small deterministic conforming manufactured meshes."""

from collections import Counter

import numpy as np

from cfd_adapt.mesh import Mesh, triangle_edges


def with_boundary(xy, triangles, name="wall"):
    edges = Counter(edge for triangle in triangles for edge in triangle_edges(triangle))
    boundary = [edge for edge, count in edges.items() if count == 1]
    return Mesh(xy, triangles, {name: boundary})


def grid(n, jitter=0.0, hole=False):
    """Unit square with n cells per side and optional central square hole."""
    if hole and n % 4:
        raise ValueError("hole requires a cell count divisible by four")
    line = np.linspace(0, 1, n + 1)
    xy = np.column_stack([v.ravel() for v in np.meshgrid(line, line)])
    interior = np.all((xy > 0) & (xy < 1), axis=1)
    if hole:
        x, y = xy.T
        inner_rim = ((np.isin(x, [0.25, 0.75]) & (y >= 0.25) & (y <= 0.75))
                     | (np.isin(y, [0.25, 0.75]) & (x >= 0.25) & (x <= 0.75)))
        interior &= ~inner_rim
    xy[interior] += np.random.default_rng(41).uniform(-jitter / n, jitter / n, (interior.sum(), 2))
    triangles = []
    for j in range(n):
        for i in range(n):
            if hole and n // 4 <= i < n * 3 // 4 and n // 4 <= j < n * 3 // 4:
                continue
            a = j * (n + 1) + i
            b, c, d = a + 1, a + n + 2, a + n + 1
            triangles.extend(((a, b, c), (a, c, d)))
    used = np.unique(triangles)
    return with_boundary(xy[used], np.searchsorted(used, triangles))
