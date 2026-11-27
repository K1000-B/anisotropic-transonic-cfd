"""Small JSON summaries and optional NumPy/ParaView diagnostic fields."""

from dataclasses import asdict
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np

from .adaptation import AdaptationResult
from .mesh import Mesh


def distribution(values):
    quantiles = np.percentile(values, [0, 5, 50, 95, 100])
    return dict(zip(("min", "p05", "median", "p95", "max"), map(float, quantiles)))


def summarise(source: Mesh, result: AdaptationResult, metric_parameters, remesh_parameters):
    """Return JSON-ready complexity, geometry, quality and operation evidence."""
    remesh = result.diagnostics
    parameters = dict(metric=asdict(metric_parameters), remeshing=asdict(remesh_parameters))

    def mesh_statistics(mesh, lengths, physical, metric, outside, penalty):
        return {
            "vertices": len(mesh.coordinates), "triangles": len(mesh.triangles),
            "edges": len(mesh.edges), "area": float(mesh.triangle_areas.sum()),
            "components": mesh.components, "boundary_components": mesh.boundary_components,
            "euler_characteristic": mesh.euler_characteristic,
            "marker_edges": {name: len(edges) for name, edges in mesh.markers.items()},
            "metric_edge_lengths": distribution(lengths), "physical_quality": distribution(physical),
            "metric_quality": distribution(metric), "short_edges": outside[0], "long_edges": outside[1],
            "mean_squared_log_length": penalty,
        }

    return {
        "parameters": parameters,
        "initial": mesh_statistics(source, remesh.initial_edge_lengths, remesh.initial_physical_quality,
                                   remesh.initial_metric_quality, remesh.initial_outside, remesh.initial_edge_penalty),
        "adapted": mesh_statistics(result.mesh, remesh.final_edge_lengths, remesh.final_physical_quality,
                                   remesh.final_metric_quality, remesh.final_outside, remesh.final_edge_penalty),
        "metric": {
            "target_complexity": float(metric_parameters.complexity),
            "optimal_complexity": result.metric.optimal_complexity,
            "bounded_complexity": result.metric.bounded_complexity,
            "complexity_relative_error": result.metric.bounded_complexity / metric_parameters.complexity - 1,
            "feasible_complexities": [value if np.isfinite(value) else None
                                      for value in result.metric.feasible_complexities],
            "log_global_scaling": result.metric.log_global_scaling, "log_beta": result.metric.log_beta,
            "bisection_iterations": result.metric.bisection_iterations,
            "anisotropy": distribution(result.metric.anisotropy),
            "target_long_size": distribution(result.metric.target_sizes[:, 0]),
            "target_short_size": distribution(result.metric.target_sizes[:, 1]),
        },
        "operations": {"accepted": dict(remesh.accepted), "rejected": dict(remesh.rejected),
                       "passes": remesh.passes, "termination": remesh.termination,
                       "constrained_noncompliant_edges": remesh.final_constrained_edges},
        "timings": result.timings,
    }


def _array(element, name, values, vtk_type="Float64"):
    attrs = {"type": vtk_type, "format": "ascii"}
    if name:
        attrs["Name"] = name
    values = np.asarray(values)
    if values.ndim > 1:
        values = values.reshape(len(values), -1)
        attrs["NumberOfComponents"] = str(values.shape[1])
    node = ET.SubElement(element, "DataArray", attrs)
    node.text = "\n" + " ".join(map(str, values.ravel())) + "\n"


def write_vtu(path, coordinates, cells, *, point_data=None, cell_data=None):
    """Write an ASCII VTK UnstructuredGrid of planar triangles or lines."""
    node_count = len(coordinates)
    n_cells, width = cells.shape
    cell_type = {2: 3, 3: 5}[width]
    root = ET.Element("VTKFile", type="UnstructuredGrid", version="0.1", byte_order="LittleEndian")
    grid = ET.SubElement(root, "UnstructuredGrid")
    piece = ET.SubElement(grid, "Piece", NumberOfPoints=str(node_count), NumberOfCells=str(n_cells))
    for name, fields, size in (("PointData", point_data, node_count), ("CellData", cell_data, n_cells)):
        if fields:
            element = ET.SubElement(piece, name)
            for label, values in fields.items():
                if len(values) != size:
                    raise ValueError(f"invalid {name} array length for {label}")
                vtk_type = "Int64" if np.asarray(values).dtype.kind in "iu" else "Float64"
                _array(element, label, values, vtk_type)
    points = ET.SubElement(piece, "Points")
    _array(points, None, np.column_stack((coordinates, np.zeros(node_count))))
    connectivity = ET.SubElement(piece, "Cells")
    _array(connectivity, "connectivity", cells.ravel(), "Int64")
    _array(connectivity, "offsets", width * np.arange(1, n_cells + 1), "Int64")
    _array(connectivity, "types", np.full(n_cells, cell_type), "UInt8")
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def _vector3(values):
    return np.pad(values, ((0, 0), (0, 1)))


def _tensor3(values):
    return np.pad(values, ((0, 0), (0, 1), (0, 1)))


def write_diagnostics(directory, source, result, summary, *, level="summary"):
    """Write summary, optional ParaView VTUs, and optional raw NumPy arrays."""
    if level not in {"summary", "vtk", "full"}:
        raise ValueError("diagnostic level must be summary, vtk or full")
    path = Path(directory)
    path.mkdir(parents=True, exist_ok=True)
    (path / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    if level == "summary":
        return
    recovery, metric, remesh = result.recovery, result.metric, result.diagnostics
    source_fields = {
        "density": result.density, "nodal_gradient": _vector3(recovery.nodal_gradient),
        "raw_hessian": _tensor3(recovery.raw_hessian), "hessian": _tensor3(recovery.hessian),
        "hessian_eigenvalues": metric.hessian_eigenvalues, "curvature": _tensor3(metric.curvature),
        "optimal_metric": _tensor3(metric.optimal), "bounded_metric": _tensor3(metric.bounded),
        "metric_density": metric.metric_density, "target_sizes": metric.target_sizes,
        "principal_long_direction": _vector3(metric.principal_directions[:, :, 0]),
        "principal_short_direction": _vector3(metric.principal_directions[:, :, 1]),
        "anisotropy": metric.anisotropy,
    }
    for label, mesh, edges, lengths, physical, quality, fields in (
        ("source", source, remesh.initial_edges, remesh.initial_edge_lengths,
         remesh.initial_physical_quality, remesh.initial_metric_quality, source_fields),
        ("adapted", result.mesh, remesh.final_edges, remesh.final_edge_lengths,
         remesh.final_physical_quality, remesh.final_metric_quality,
         {"source_vertex_index": result.source_vertex_indices}),
    ):
        cell_fields = {"physical_quality": physical, "metric_quality": quality}
        if label == "source":
            cell_fields["element_gradient"] = _vector3(recovery.element_gradient)
        write_vtu(path / f"{label}.vtu", mesh.coordinates, mesh.triangles,
                  point_data=fields, cell_data=cell_fields)
        write_vtu(path / f"{label}_edges.vtu", mesh.coordinates, edges,
                  cell_data={"metric_edge_length": lengths})
    if level == "full":
        np.savez_compressed(
            path / "fields.npz", source_coordinates=source.coordinates,
            source_triangles=source.triangles, source_point_ids=source.point_ids,
            adapted_coordinates=result.mesh.coordinates,
            adapted_triangles=result.mesh.triangles, source_vertex_indices=result.source_vertex_indices,
            density=result.density, **vars(recovery), **{
                name: values for name, values in vars(metric).items() if isinstance(values, np.ndarray)
            }, metric_density=metric.metric_density, target_sizes=metric.target_sizes,
            anisotropy=metric.anisotropy, **{
                name: values for name, values in vars(remesh).items() if isinstance(values, np.ndarray)
            },
        )
