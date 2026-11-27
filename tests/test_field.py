import unittest

import numpy as np
from numpy.testing import assert_allclose

from cfd_adapt.field import OutsideSourceError, P1MetricField
from support import grid


class FieldTests(unittest.TestCase):
    def test_constant_and_location(self):
        mesh = grid(8)
        matrix = np.array([[4, 1], [1, 3]])
        nodal = np.broadcast_to(matrix, (len(mesh.coordinates), 2, 2)).copy()
        field = P1MetricField(mesh, nodal)
        nodal[:] = 0
        for a, b in (([0, 0], [1, 1]), ([1, 0], [0, 1]), ([0, 0.2], [0, 0.95]),
                     ([0.1, 0.6], [0.9, 0.2])):
            d = np.subtract(b, a)
            length, mid = field.edge_measure(a, b, midpoint=True)
            self.assertAlmostEqual(length, np.sqrt(d @ matrix @ d), places=13)
            self.assertEqual(length, field.edge_measure(b, a))
            assert_allclose(field.at(mid), matrix)
            assert_allclose(mid, np.add(a, b) / 2, atol=1e-12)
        assert_allclose(field.at([-1e-12, 0.5]), matrix)
        with self.assertRaises(OutsideSourceError):
            field.at([-1e-5, 0.5])

    def test_piecewise_affine_integral(self):
        mesh = grid(8)
        x, y = mesh.coordinates.T
        scale = 1 + 7 * abs(x - 0.5) + 0.2 * y
        field = P1MetricField(mesh, scale[:, None, None] * np.eye(2))
        y0 = 0.37
        a, b = np.array([0.17, y0]), np.array([0.93, y0])
        offset = 1 + 0.2 * y0
        at_centre = offset ** 1.5
        expected = 2/21 * ((offset + 7 * 0.33) ** 1.5 + (offset + 7 * 0.43) ** 1.5 - 2 * at_centre)
        length, mid = field.edge_measure(a, b, midpoint=True)
        self.assertAlmostEqual(length, expected, places=13)
        self.assertAlmostEqual(field.edge_measure(a, mid), 0.5 * expected, places=12)
        self.assertAlmostEqual(field.edge_measure(mid, b), 0.5 * expected, places=12)
        assert_allclose(field.at([0.26, y0]), (offset + 7 * 0.24) * np.eye(2))
        lengths = field.source_edge_lengths(list(mesh.edges))
        assert_allclose(lengths, [field.edge_measure(*mesh.coordinates[list(e)]) for e in mesh.edges])

    def test_hole_and_boundary(self):
        mesh = grid(8, hole=True)
        field = P1MetricField(mesh, np.tile(np.eye(2), (len(mesh.coordinates), 1, 1)))
        with self.assertRaises(OutsideSourceError):
            field.edge_measure([0.1, 0.5], [0.9, 0.5])
        self.assertAlmostEqual(field.edge_measure([0.1, 0.25], [0.9, 0.25]), 0.8)
