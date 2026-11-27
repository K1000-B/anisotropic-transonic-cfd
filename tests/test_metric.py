import unittest

import numpy as np
from numpy.testing import assert_allclose

from cfd_adapt.mesh import Mesh
from cfd_adapt.metric import (
    MetricParameters, bounded_log_spectrum, construct_metric, nodal_average,
    recover_hessian, regularise_curvature, tensor,
)
from support import grid


class MeshTests(unittest.TestCase):
    def test_mesh_and_lumped_areas(self):
        mesh = grid(8, jitter=0.1, hole=True)
        self.assertEqual((mesh.components, mesh.boundary_components, mesh.euler_characteristic), (1, 2, 0))
        self.assertAlmostEqual(mesh.nodal_areas.sum(), 0.75)
        self.assertAlmostEqual(mesh.nodal_areas.sum(), mesh.triangle_areas.sum())
        with self.assertRaises(ValueError):
            mesh.coordinates.setflags(write=True)
        invalid = mesh.triangles.copy()
        invalid[0, :2] = invalid[0, 1::-1]
        with self.assertRaisesRegex(ValueError, "counterclockwise"):
            Mesh(mesh.coordinates, invalid, mesh.markers)
        with self.assertRaisesRegex(ValueError, "markers"):
            Mesh(mesh.coordinates, mesh.triangles, {})
        invalid[0, 0] = len(mesh.coordinates)
        with self.assertRaisesRegex(ValueError, "connectivity"):
            Mesh(mesh.coordinates, invalid, mesh.markers)

    def test_area_weighting(self):
        mesh = Mesh([[0, 0], [1, 0], [1, 1], [0, 2]], [[0, 1, 2], [0, 2, 3]],
                    {"boundary": [[0, 1], [1, 2], [2, 3], [3, 0]]})
        assert_allclose(mesh.nodal_areas, [0.5, 1/6, 0.5, 1/3])
        assert_allclose(nodal_average(mesh, [1, 4]), [3, 1, 3, 4])


class RecoveryTests(unittest.TestCase):
    def test_constant_and_affine(self):
        mesh = grid(9, jitter=0.22)
        for bx, by in [(0, 0), (2, -3)]:
            result = recover_hessian(mesh, 4 + mesh.coordinates @ [bx, by])
            assert_allclose(result.element_gradient, np.tile([bx, by], (len(mesh.triangles), 1)), atol=3e-14)
            assert_allclose(result.nodal_gradient, np.tile([bx, by], (len(mesh.coordinates), 1)), atol=3e-14)
            assert_allclose(result.hessian, 0, atol=3e-13)
        with self.assertRaisesRegex(ValueError, "finite"):
            recover_hessian(mesh, np.full(len(mesh.coordinates), np.nan))

    def test_quadratic_recovery(self):
        exact = np.array([[1., 0.4], [0.4, -0.6]])
        for jitter in (0, 0.17):
            errors = []
            for n in (12, 24, 48):
                mesh = grid(n, jitter)
                rho = 1 + 0.5 * np.einsum("ni,ij,nj->n", mesh.coordinates, exact, mesh.coordinates)
                result = recover_hessian(mesh, rho)
                assert_allclose(result.hessian, result.hessian.transpose(0, 2, 1), atol=0, rtol=0)
                e = np.linalg.norm(result.hessian - exact, axis=(1, 2))
                inward = np.all((mesh.coordinates > 2.1/n) & (mesh.coordinates < 1 - 2.1/n), axis=1)
                rim = np.unique(list(mesh.boundary))
                self.assertLess(e[inward].mean(), 0.12 if jitter else 1e-10)
                self.assertGreater(e[rim].mean(), 0.28)
                errors.append(np.sqrt(np.average(e * e, weights=mesh.nodal_areas)))
            self.assertLess(errors[-1], 0.8 * errors[0])

    def test_spectral_floor(self):
        h = np.array([[5, 2], [-5, -2], [5, -2], [5, 0], [0, 0], [-3, -3]])
        q = np.array([[0.6, -0.8], [0.8, 0.6]])
        matrices = np.stack([q @ np.diag(values) @ q.T for values in h])
        matrices = (matrices + matrices.transpose(0, 2, 1)) / 2
        signed, basis, floored = regularise_curvature(matrices, 0.5)
        expected = np.stack([q @ np.diag(np.maximum(abs(values), 0.5)) @ q.T for values in h])
        assert_allclose(tensor(basis, signed), matrices, atol=1e-14)
        assert_allclose(tensor(basis, floored), expected, atol=1e-14)


class MetricTests(unittest.TestCase):
    def test_optimal_and_bounded(self):
        mesh = grid(4)
        q = np.array([[0.6, -0.8], [0.8, 0.6]])
        eigenvalues = np.column_stack((np.full(len(mesh.coordinates), 0.2),
                                       2 + 70 * mesh.coordinates[:, 0]))
        h = tensor(np.broadcast_to(q, (len(eigenvalues), 2, 2)), eigenvalues)
        settings = MetricParameters(0.1, 0.05, 0.4, 4, p=2, target_complexity=45)
        result = construct_metric(mesh, h, settings)
        det = np.prod(eigenvalues, axis=1)
        s = 45 / np.sum(mesh.nodal_areas * det ** (1/3))
        assert_allclose(result.optimal, (s * det ** (-1/6))[:, None, None] * h)
        self.assertAlmostEqual(result.optimal_complexity, 45)
        self.assertGreater(result.bisection_iterations, 0)
        self.assertNotEqual(result.log_beta, 0)
        self.assertAlmostEqual(result.bounded_complexity, 45, places=7)
        matrices = np.einsum("ia,nij,jb->nab", q, result.bounded, q)
        assert_allclose(matrices[:, 0, 1], 0, atol=1e-12)
        numerical_complexity = mesh.nodal_areas @ np.sqrt(np.linalg.det(result.bounded))
        self.assertAlmostEqual(numerical_complexity, 45, places=7)
        self.assertGreaterEqual(result.target_sizes.min(), 0.05 - 1e-14)
        self.assertLessEqual(result.target_sizes.max(), 0.4 + 1e-14)
        self.assertLessEqual(result.anisotropy.max(), 4 + 1e-14)

    def test_exact_order_and_bounded_plateau(self):
        logs = bounded_log_spectrum(np.log([[0.001, 500]]), 0, np.log(1), np.log(64), np.log(2))
        assert_allclose(np.exp(logs), [[16, 64]])
        mesh = grid(2)
        h = np.zeros((len(mesh.coordinates), 2, 2))
        h[0] = np.diag([1e-110, 1e110])
        for target, value in ((4, 4), (100, 100)):
            settings = MetricParameters(1e-100, 0.1, 0.5, 4, target_complexity=target)
            result = construct_metric(mesh, h, settings)
            assert_allclose(result.bounded_eigenvalues, value)
        with self.assertRaisesRegex(ValueError, "infeasible"):
            construct_metric(mesh, h, MetricParameters(1e-100, 0.1, 0.5, 4, target_complexity=3))
        with self.assertRaisesRegex(ValueError, "positive integer"):
            MetricParameters(0.1, 0.1, 1, 1, p=1.5, target_nodes=2)
