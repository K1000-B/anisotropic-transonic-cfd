"""SU2 ASCII mesh and nodal CSV interfaces for single-zone planar triangles."""

import csv
from pathlib import Path

import numpy as np

from .mesh import Mesh, cross
from .remesh import GeometryConstraints


def read_su2(path):
    """Read planar SU2 type-5 triangles and type-3 marker segments.

    Explicit point IDs are retained and must be present for either every
    coordinate row or none. Clockwise elements are oriented counterclockwise.
    """
    with Path(path).open(encoding="utf-8-sig") as stream:
        lines = iter(line for raw in stream if (line := raw.split("%", 1)[0].strip()))
        sections = {}
        while (line := next(lines, None)) is not None:
            key, sep, value = line.partition("=")
            key, value = key.strip(), value.strip()
            if not sep or key not in {"NDIME", "NELEM", "NPOIN", "NMARK"} or key in sections:
                raise ValueError(f"invalid or unsupported SU2 section {line!r}")
            parts = value.split()
            if len(parts) != 1:
                raise ValueError(f"unsupported SU2 count {line!r}")
            count = int(parts[0])
            if count < 0:
                raise ValueError(f"negative SU2 count in {line!r}")
            if key == "NDIME":
                if count != 2:
                    raise ValueError("only planar SU2 meshes (NDIME=2) are supported")
                sections[key] = count
            elif key == "NMARK":
                markers = {}
                for _ in range(count):
                    name = _assignment(lines, "MARKER_TAG")
                    marker_count = int(_assignment(lines, "MARKER_ELEMS"))
                    if name in markers or marker_count < 1:
                        raise ValueError("duplicate or empty SU2 marker")
                    markers[name] = [_element(lines, "3", 2, False) for _ in range(marker_count)]
                sections[key] = markers
            elif key == "NELEM":
                sections[key] = [_element(lines, "5", 3, True) for _ in range(count)]
            else:
                sections[key] = [_row(lines) for _ in range(count)]
    if set(sections) != {"NDIME", "NELEM", "NPOIN", "NMARK"}:
        raise ValueError("missing required SU2 sections")
    nodes = sections["NPOIN"]
    widths = {len(row) for row in nodes}
    if widths not in ({2}, {3}) or not sections["NELEM"]:
        raise ValueError("inconsistent point rows or empty SU2 mesh")
    coordinates = np.array([[float(row[0]), float(row[1])] for row in nodes])
    point_ids = [int(row[2]) for row in nodes] if widths == {3} else list(range(len(nodes)))
    index = {point_id: i for i, point_id in enumerate(point_ids)}
    if len(index) != len(nodes):
        raise ValueError("duplicate SU2 point IDs")
    try:
        triangles = np.array([[index[v] for v in tri] for tri in sections["NELEM"]])
        markers = {name: [[index[v] for v in edge] for edge in segments]
                   for name, segments in sections["NMARK"].items()}
    except KeyError as error:
        raise ValueError(f"unknown SU2 point ID in connectivity {error}") from error
    a, b, c = coordinates[triangles].transpose(1, 0, 2)
    clockwise = cross(b - a, c - a) < 0
    triangles[clockwise] = triangles[clockwise][:, [0, 2, 1]]
    return Mesh(coordinates, triangles, markers, np.array(point_ids))


def _assignment(lines, expected):
    line = next(lines, "")
    key, sep, value = line.partition("=")
    if not sep or key.strip() != expected or not value.strip():
        raise ValueError(f"expected {expected}, got {line!r}")
    return value.strip()


def _row(lines):
    row = next(lines, "").split()
    if not row or "=" in row:
        raise ValueError("truncated or malformed SU2 section")
    return row


def _element(lines, code, width, optional_id):
    row = _row(lines)
    if row[0] != code or len(row) not in (width + 1, width + 1 + optional_id):
        raise ValueError(f"expected SU2 element type {code} with {width} vertices")
    # Parse an optional cell ID too, although volume order defines output IDs.
    return [int(value) for value in row[1:]][:width]


def write_su2(path, mesh: Mesh):
    """Write SU2 ASCII with contiguous output point and triangle IDs."""
    for name in mesh.markers:
        if "%" in name or "\r" in name:
            raise ValueError(f"marker {name!r} cannot be represented in SU2 ASCII")
    with Path(path).open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(f"NDIME= 2\nNELEM= {len(mesh.triangles)}\n")
        for i, (a, b, c) in enumerate(mesh.triangles):
            stream.write(f"5 {a} {b} {c} {i}\n")
        stream.write(f"NPOIN= {len(mesh.coordinates)}\n")
        for i, (x, y) in enumerate(mesh.coordinates):
            stream.write(f"{x:.17g} {y:.17g} {i}\n")
        stream.write(f"NMARK= {len(mesh.markers)}\n")
        for name, segments in mesh.markers.items():
            stream.write(f"MARKER_TAG= {name}\nMARKER_ELEMS= {len(segments)}\n")
            for a, b in segments:
                stream.write(f"3 {a} {b}\n")


def read_density_csv(path, mesh: Mesh, *, density_column="Density", id_column="PointID",
                     coordinate_columns=("x", "y"), coordinate_atol=None):
    """Align a shuffled nodal CSV to mesh vertices using IDs and/or coordinates.

    With both kinds of identity, validate coordinates against each matched ID.
    With coordinates alone, each row must match exactly one vertex within the
    absolute componentwise tolerance (default 1e-9 times the largest extent).
    """
    if coordinate_atol is None:
        coordinate_atol = 1e-9 * float(np.ptp(mesh.coordinates, axis=0).max())
    if not np.isfinite(coordinate_atol) or coordinate_atol < 0:
        raise ValueError("coordinate tolerance must be finite and nonnegative")
    output = np.empty(len(mesh.coordinates))
    by_id = {int(point_id): i for i, point_id in enumerate(mesh.point_ids)}
    order = np.argsort(mesh.coordinates[:, 0])
    x = mesh.coordinates[order, 0]
    seen = set()
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, skipinitialspace=True)
        header = [name.strip() for name in next(reader, [])]
        if not header or len(header) != len(set(header)) or density_column not in header:
            raise ValueError(f"CSV requires a unique {density_column!r} column")
        column = {name: i for i, name in enumerate(header)}
        has_id = id_column in column
        has_xy = all(name in column for name in coordinate_columns)
        if (not has_id and not has_xy) or (not has_xy and any(name in column for name in coordinate_columns)):
            raise ValueError("CSV needs PointID, x and y, or all three")
        for row in reader:
            if len(row) != len(header):
                raise ValueError(f"CSV row {reader.line_num} has the wrong number of fields")
            density = float(row[column[density_column]])
            if not np.isfinite(density):
                raise ValueError(f"nonfinite density on CSV row {reader.line_num}")
            if has_xy:
                coords = np.array([float(row[column[name]]) for name in coordinate_columns])
                if not np.isfinite(coords).all():
                    raise ValueError(f"nonfinite coordinates on CSV row {reader.line_num}")
            if has_id:
                point_id = int(row[column[id_column]])
                if point_id not in by_id:
                    raise ValueError(f"unknown solution point ID {point_id}")
                index = by_id[point_id]
                if has_xy and np.any(np.abs(coords - mesh.coordinates[index]) > coordinate_atol):
                    raise ValueError(f"solution coordinates disagree with mesh point ID {point_id}")
            else:
                lo = np.searchsorted(x, coords[0] - coordinate_atol, side="left")
                hi = np.searchsorted(x, coords[0] + coordinate_atol, side="right")
                candidates = order[lo:hi]
                match = candidates[np.all(np.abs(mesh.coordinates[candidates] - coords) <= coordinate_atol, axis=1)]
                if len(match) != 1:
                    raise ValueError(f"CSV row {reader.line_num} has {len(match)} coordinate matches")
                index = int(match[0])
            if index in seen:
                raise ValueError(f"multiple solution rows match mesh vertex {index}")
            seen.add(index)
            output[index] = density
    if len(seen) != len(output):
        raise ValueError(f"solution is missing {len(output) - len(seen)} mesh vertices")
    return output


def constraints_from_ids(mesh: Mesh, *, fixed_point_ids=(), protected_point_id_edges=(),
                         splittable_point_id_edges=()):
    """Translate explicit external IDs from the loaded SU2 mesh to array indices."""
    index = {int(point_id): i for i, point_id in enumerate(mesh.point_ids)}

    def edges(data):
        return frozenset((index[a], index[b]) for a, b in data)

    try:
        return GeometryConstraints(
            frozenset(index[v] for v in fixed_point_ids), edges(protected_point_id_edges),
            None if splittable_point_id_edges is None else edges(splittable_point_id_edges),
        )
    except KeyError as error:
        raise ValueError(f"unknown constraint point ID {error}") from error
