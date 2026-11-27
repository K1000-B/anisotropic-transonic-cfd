import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from cfd_adapt import GeometryConstraints, MetricParameters, RemeshParameters, adapt
from cfd_adapt.field import P1MetricField
from cfd_adapt.mesh import edge_key
from cfd_adapt.remesh import _Triangulation, remesh, triangle_quality
from cfd_adapt.verification import grid, with_boundary


def state(mesh, scale=1, constraints=GeometryConstraints()):
    return _Triangulation(P1MetricField(mesh, np.tile(scale * np.eye(2),
                                                      (len(mesh.coordinates), 1, 1))), constraints)


class OperationTests(unittest.TestCase):
    def test_metric_quality_equilateral(self):
        xy = np.array([[[0, 0], [1, 0], [0.5, 0.5 * np.sqrt(3)]]])
        assert_allclose(triangle_quality(xy), 1)
        scaled = xy * np.array([0.5, 2])
        assert_allclose(triangle_quality(scaled, np.array([np.diag([4, 0.25])])), 1)

    def test_interior_boundary_split_and_markers(self):
        mesh = grid(1)
        split = state(mesh, 4, GeometryConstraints(splittable_boundary_edges=frozenset({(0, 1)})))
        self.assertEqual(split.split((0, 2)), "constrained")
        self.assertIsNone(split.split((0, 3)))
        assert_allclose(split.xy[4], [0.5, 0.5], atol=1e-12)
        self.assertIsNone(split.split((0, 1)))
        result, _ = split.snapshot()
        self.assertEqual((len(result.coordinates), len(result.triangles)), (6, 5))
        self.assertEqual(result.boundary[0, 5], "wall")
        self.assertEqual(result.boundary[1, 5], "wall")
        self.assertEqual(result.euler_characteristic, 1)
        self.assertAlmostEqual(result.triangle_areas.sum(), 1)

    def test_protected_collapse_and_local_coarsening(self):
        mesh = grid(4)
        edge = (6, 7)
        config = RemeshParameters()
        fixed = state(mesh, constraints=GeometryConstraints(protected_edges=frozenset({edge})))
        self.assertEqual(fixed.collapse(edge, config), "constrained")
        self.assertEqual(fixed.relocate(6, config), "constrained")
        self.assertEqual(fixed.swap(edge, config), "constrained")
        self.assertEqual(fixed.split(edge), "constrained")
        editable = state(mesh)
        self.assertEqual(editable.collapse((0, 1), config), "constrained")
        self.assertIsNone(editable.collapse(edge, config))
        result, origin = editable.snapshot()
        self.assertEqual(len(result.coordinates), len(mesh.coordinates) - 1)
        self.assertAlmostEqual(result.triangle_areas.sum(), 1)
        for index, original in enumerate(origin):
            assert_array_equal(result.coordinates[index], mesh.coordinates[original])

    def test_swap_improves_quality_and_rejects_concave_patch(self):
        mesh = with_boundary([[0, 0], [1, 0], [1, 1], [0, 0.25]], [[0, 1, 2], [0, 2, 3]])
        editable = state(mesh)
        before = editable.quality(list(editable.triangles.values())).min()
        self.assertIsNone(editable.swap((0, 2), RemeshParameters()))
        self.assertIn((1, 3), editable.edges)
        self.assertGreater(editable.quality(list(editable.triangles.values())).min(), before)
        self.assertEqual(editable.swap((1, 3), RemeshParameters()), "quality")
        concave = with_boundary([[0, 0], [1, 0], [0.25, 0.25], [0, 1]], [[0, 1, 2], [0, 2, 3]])
        self.assertEqual(state(concave).swap((0, 2), RemeshParameters()), "orientation")

    def test_relocation_and_inverted_trial_rejection(self):
        triangles = [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]
        mesh = with_boundary([[0, 0], [1, 0], [1, 1], [0, 1], [0.25, 0.4]], triangles)
        editable = state(mesh)
        config = RemeshParameters(collapse_below=0.3)
        before = editable.quality(list(editable.triangles.values())).min()
        self.assertIsNone(editable.relocate(4, config))
        assert_allclose(editable.xy[4], [0.5, 0.5])
        self.assertGreater(editable.quality(list(editable.triangles.values())).min(), before)
        self.assertEqual(editable.relocate(0, config), "constrained")
        # This centre sees a re-entrant boundary. Its neighbour mean lies
        # beyond the admissible patch, so the one-trial search must reject it.
        mesh = with_boundary([[0, 0], [2, 0], [2, 2], [1, 0.2], [0, 2], [1, 0.1]],
                             [[0, 1, 5], [1, 2, 5], [2, 3, 5], [3, 4, 5], [4, 0, 5]])
        editable = state(mesh)
        rejection = editable.relocate(5, RemeshParameters(relocation_steps=1))
        self.assertIn(rejection, {"orientation", "outside_source"})
        assert_array_equal(editable.xy[5], [1, 0.1])

    def test_hole_and_operation_limit(self):
        mesh = grid(8, jitter=0.2, hole=True)
        field = P1MetricField(mesh, np.tile(100 * np.eye(2), (len(mesh.coordinates), 1, 1)))
        output, diag, _ = remesh(field, RemeshParameters(max_operations=3),
                                 GeometryConstraints(splittable_boundary_edges=None))
        self.assertEqual(diag.termination, "operation_limit")
        self.assertEqual(diag.accepted.total(), 3)
        self.assertEqual((output.euler_characteristic, output.boundary_components), (0, 2))
        self.assertAlmostEqual(output.triangle_areas.sum(), 0.75)

    def test_stationary_unmet_and_criteria_met(self):
        mesh = with_boundary([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]], [[0, 1, 2]])
        for scale, expected in ((1, "criteria_met"), (4, "stationary")):
            field = P1MetricField(mesh, np.tile(scale * np.eye(2), (3, 1, 1)))
            _, diag, _ = remesh(field)
            self.assertEqual(diag.termination, expected)
            self.assertEqual(diag.final_constrained_edges, 3 if scale == 4 else 0)


class IntegrationTests(unittest.TestCase):
    def test_rotated_manufactured_fields(self):
        mesh = grid(12, jitter=0.15)
        for angle in (19, 67):
            theta = np.deg2rad(angle)
            n = np.array([np.cos(theta), np.sin(theta)])
            t = np.array([-n[1], n[0]])
            x = (mesh.coordinates - 0.5) @ n
            density = 1 + 0.3 * (1 + np.tanh(x / 0.09))
            result = adapt(mesh, density, MetricParameters(0.05, 0.023, 0.35, 7, target_nodes=180),
                           RemeshParameters(max_passes=5),
                           GeometryConstraints(splittable_boundary_edges=None))
            self.assertAlmostEqual(result.metric.bounded_complexity, np.sqrt(3) * 90, places=6)
            band = (abs(x) >= 0.045) & (abs(x) <= 0.11)
            bounded = result.metric.bounded[band]
            self.assertGreater(np.median((n @ bounded @ n) / (t @ bounded @ t)), 8)
            self.assertGreater(result.diagnostics.accepted["split"], 0)
            self.assertGreater(result.diagnostics.accepted["collapse"], 0)
            self.assertLess(result.diagnostics.final_edge_penalty, result.diagnostics.initial_edge_penalty)
            self.assertGreater(np.median(result.diagnostics.final_metric_quality), 0.75)
            self.assertAlmostEqual(result.mesh.triangle_areas.sum(), 1)
            self.assertEqual(result.mesh.euler_characteristic, 1)
            assert_array_equal(result.density, density)
