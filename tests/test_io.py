import csv
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
import xml.etree.ElementTree as ET

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

from cfd_adapt import Mesh, MetricParameters, RemeshParameters, adapt
from cfd_adapt.diagnostics import summarise, write_diagnostics
from cfd_adapt.io import constraints_from_ids, read_density_csv, read_su2, write_su2
from support import grid


class IOTests(unittest.TestCase):
    def test_mission_01_round_trip(self):
        baseline = Path(__file__).resolve().parents[1] / "01_baseline_mesh/reference/baseline.su2"
        source = read_su2(baseline)
        self.assertEqual(len(source.coordinates), 7154)
        self.assertEqual(len(source.triangles), 13702)
        self.assertEqual(source.euler_characteristic, 0)
        self.assertEqual({name: len(edges) for name, edges in source.markers.items()},
                         {"airfoil": 524, "farfield": 62, "outlet": 20})
        with TemporaryDirectory() as directory:
            export = Path(directory) / "mesh.su2"
            write_su2(export, source)
            reread = read_su2(export)
        assert_array_equal(source.coordinates, reread.coordinates)
        assert_array_equal(source.triangles, reread.triangles)
        self.assertEqual(source.boundary, reread.boundary)

    def test_csv_association(self):
        mesh = grid(3)
        ids = 8 + 2 * np.arange(len(mesh.coordinates))
        mesh = Mesh(mesh.coordinates, mesh.triangles, mesh.markers, ids)
        density = 1 + mesh.coordinates @ [2, 3]
        rows = [[i, *xy, value] for i, xy, value in zip(ids, mesh.coordinates, density)][::-1]
        with TemporaryDirectory() as directory:
            path = Path(directory) / "solution.csv"

            def write(header, rows):
                with path.open("w", newline="") as stream:
                    csv.writer(stream).writerows([header, *rows])

            header = ["PointID", "x", "y", "Density"]
            for columns in ([0, 1, 2, 3], [0, 3], [1, 2, 3]):
                write([header[i] for i in columns], [[row[i] for i in columns] for row in rows])
                assert_array_equal(read_density_csv(path, mesh), density)
            write(header, rows[:-1])
            with self.assertRaisesRegex(ValueError, "missing"):
                read_density_csv(path, mesh)
            write(header, rows + rows[:1])
            with self.assertRaisesRegex(ValueError, "multiple"):
                read_density_csv(path, mesh)
            write(header, [[rows[0][0], *rows[1][1:3], rows[0][3]], *rows[1:]])
            with self.assertRaisesRegex(ValueError, "disagree"):
                read_density_csv(path, mesh)
            write(header[1:], [row[1:] for row in rows])
            with self.assertRaisesRegex(ValueError, "coordinate matches"):
                read_density_csv(path, mesh, coordinate_atol=0.4)

    def test_explicit_mesh_point_ids(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "explicit.su2"
            path.write_text("NDIME=2\nNELEM=1\n5 6 9 4 0\nNPOIN=3\n0 0 6\n1 0 4\n0 1 9\n"
                            "NMARK=1\nMARKER_TAG=wall\nMARKER_ELEMS=3\n3 6 4\n3 4 9\n3 9 6\n")
            mesh = read_su2(path)
        assert_array_equal(mesh.point_ids, [6, 4, 9])
        assert_array_equal(mesh.triangles, [[0, 1, 2]])
        constraints = constraints_from_ids(mesh, fixed_point_ids=[6], splittable_point_id_edges=[(4, 9)])
        self.assertEqual(constraints.fixed_vertices, {0})
        self.assertEqual(constraints.splittable_boundary_edges, {(1, 2)})

    def test_diagnostics(self):
        mesh = grid(4)
        settings = MetricParameters(0.1, 0.05, 0.8, 3, target_nodes=25)
        passes = RemeshParameters(max_passes=1)
        result = adapt(mesh, mesh.coordinates[:, 0] ** 2, settings, passes)
        summary = summarise(mesh, result, settings, passes)
        with TemporaryDirectory() as directory:
            path = Path(directory)
            write_diagnostics(path, mesh, result, summary)
            self.assertEqual([x.name for x in path.iterdir()], ["summary.json"])
            write_diagnostics(path, mesh, result, summary, level="full")
            self.assertEqual(ET.parse(path / "source.vtu").getroot().attrib["type"], "UnstructuredGrid")
            with np.load(path / "fields.npz") as fields:
                self.assertEqual(fields["bounded"].shape, (len(mesh.coordinates), 2, 2))
                assert_allclose(fields["density"], result.density)
                self.assertEqual({tuple(edge) for edge in fields["final_edges"]}, set(result.mesh.edges))
