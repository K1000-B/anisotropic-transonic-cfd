"""
SU2 Mesh and Scalar Field Extraction
====================================

This program extracts a two-dimensional triangular mesh and a
nodal scalar field from SU2 simulation files.

INPUTS
------
1. SU2 ASCII mesh file (.su2)
2. SU2 ASCII restart file (.csv)
3. Name of the scalar field to extract
4. Output file path (.npz)

OUTPUT
------
A NumPy NPZ file containing seven arrays:

    points          : (N_vertices, 2)
    triangles       : (N_triangles, 3)
    scalar          : (N_vertices,)
    boundary_edges  : (N_boundary_edges, 2)
    boundary_tags   : (N_boundary_edges,)
    boundary_nodes  : (N_boundary_nodes,)
    field_name      : scalar string

FUNDAMENTAL CONVENTION
----------------------
The scalar value stored at scalar[i] corresponds to the physical
vertex whose coordinates are stored at points[i].

All connectivity indices directly reference the points array.
Boundary nodes identify the fixed vertices on all physical boundaries.

Consequently:

    points[triangles] -> coordinates of all triangle vertices
    scalar[triangles] -> scalar values at all triangle vertices

This program does not perform any mesh adaptation, interpolation,
physical computation or modification of the original mesh.

DEPENDENCIES
------------
- NumPy
- Python standard library

USAGE
-----
python extract_solution.py mesh.su2 restart_flow.csv Mach input.npz
"""

import argparse
import csv
from pathlib import Path

import numpy as np

# ============================================================
# 1. READ THE SU2 MESH
# ============================================================

def read_su2_mesh(mesh_path):
    """
    Read a two-dimensional triangular SU2 ASCII mesh.

    The SU2 mesh format is divided into several sections:

        NDIME : spatial dimension
        NELEM : volume element connectivity
        NPOIN : vertex coordinates
        NMARK : boundary marker definitions

    SU2 uses VTK element-type identifiers:

        3 : linear edge
        5 : linear triangle
        9 : quadrilateral

    Only linear triangles and linear boundary edges are supported.

    Parameters
    ----------
    mesh_path : str or Path
        Path to the input SU2 mesh file.

    Returns
    -------
    points : ndarray of float64, shape (N_vertices, 2)
        Cartesian coordinates of the mesh vertices.

    triangles : ndarray of int64, shape (N_triangles, 3)
        Indices of the three vertices of every triangle.

    boundary_edges : ndarray of int64, shape (N_boundary_edges, 2)
        Indices of the two endpoints of every boundary edge.

    boundary_tags : ndarray of str, shape (N_boundary_edges,)
        Physical boundary marker associated with each boundary edge.

    Notes
    -----
    The original vertex numbering, element ordering, coordinates
    and triangle orientations are preserved.
    """

    # Initialise the mesh variables.
    # They will be filled when their corresponding SU2 sections
    # are encountered during file parsing.

    dimension = None
    points = None
    triangles = None

    # Boundary edges are distributed across several marker blocks.
    # Lists allow us to collect them progressively before converting
    # them into NumPy arrays.

    boundary_edges = []
    boundary_tags = []

    # Store the names of the sections already encountered.
    # This prevents accidental duplicate section definitions.

    sections = set()

    # --------------------------------------------------------
    # 1.1. Read and clean the input file
    # --------------------------------------------------------

    # SU2 uses the "%" character to introduce comments.
    #
    # Example:
    #
    #     NDIME= 2  % Spatial dimension
    #
    # Only the part preceding "%" is relevant to the parser.
    #
    # strip() also removes leading/trailing whitespace and
    # newline characters.

    with open(mesh_path, "r", encoding="utf-8") as file:

        lines = [
            line.split("%", 1)[0].strip()
            for line in file
        ]

    # Remove empty lines and create an iterator.
    #
    # An iterator allows us to consume records sequentially,
    # without manually maintaining a global line index.

    lines = iter(line for line in lines if line)

    # --------------------------------------------------------
    # 1.2. Local parsing utilities
    # --------------------------------------------------------

    def read_line():
        """
        Return the next nonempty mesh record.

        An unexpected end of the file indicates that a section
        contains fewer records than declared.
        """

        line = next(lines, None)

        if line is None:
            raise ValueError(
                "Unexpected end of SU2 mesh file."
            )

        return line

    def read_keyword(line, expected=None):
        """
        Split a SU2 section header into its keyword and value.

        Example
        -------
        Input:

            "NELEM= 4000"

        Output:

            ("NELEM", "4000")

        The optional expected argument allows us to check that
        a specific keyword appears at the appropriate location.
        """

        key, separator, value = line.partition("=")

        # Every section header must contain an equal sign.

        if not separator:
            raise ValueError(
                f"Invalid mesh record: {line}"
            )

        key = key.strip()
        value = value.strip()

        # Some sections contain mandatory nested keywords.
        # For instance, MARKER_TAG must be followed by
        # MARKER_ELEMS.

        if expected is not None and key != expected:
            raise ValueError(
                f"Expected {expected}, found {key}."
            )

        return key, value

    def read_count(value):
        """
        Convert a section record count to a nonnegative integer.

        Examples:

            NPOIN= 5000
            NELEM= 9700
            MARKER_ELEMS= 300
        """

        count = int(value)

        if count < 0:
            raise ValueError("Negative record count.")

        return count

    # --------------------------------------------------------
    # 1.3. Parse the mesh sections
    # --------------------------------------------------------

    # Read the file sequentially.
    #
    # Each section is identified by its keyword rather than
    # its absolute position in the file.
    #
    # Therefore, NELEM may appear before or after NPOIN.

    for line in lines:

        key, value = read_keyword(line)

        # A complete unpartitioned mesh must not contain
        # multiple definitions of the same main section.

        if key in sections:
            raise ValueError(
                f"Duplicate mesh section: {key}"
            )

        sections.add(key)

        # ====================================================
        # NDIME — SPATIAL DIMENSION
        # ====================================================

        # Example:
        #
        #     NDIME= 2
        #
        # This specifies that every vertex has two physical
        # coordinates, x and y.

        if key == "NDIME":

            dimension = int(value)

            # The future adaptation algorithm is formulated
            # exclusively for two-dimensional meshes.

            if dimension != 2:
                raise ValueError(
                    "Only 2D meshes are supported."
                )

        # ====================================================
        # NELEM — VOLUME ELEMENT CONNECTIVITY
        # ====================================================

        # Example:
        #
        #     NELEM= 2
        #     5 0 1 2 0
        #     5 0 2 3 1
        #
        # Each triangular element record contains:
        #
        #     type vertex_1 vertex_2 vertex_3 [element_id]
        #
        # The first integer is the element type.
        # The next three integers are the vertex indices.
        # The optional final integer is the element index.
        #
        # Only the three vertex indices are required by our
        # mesh adaptation algorithm.

        elif key == "NELEM":

            n_triangles = read_count(value)

            # Preallocate the connectivity array.
            #
            # Every row corresponds to one triangle.
            # Every column contains one of its vertex indices.
            #
            # For triangle T:
            #
            #     triangles[T] = [i, j, l]
            #
            # Its geometric coordinates can later be obtained
            # directly using:
            #
            #     points[triangles[T]]

            triangles = np.empty(
                (n_triangles, 3),
                dtype=np.int64
            )

            for i in range(n_triangles):

                # Split the element record into individual
                # whitespace-separated values.

                values = read_line().split()

                # A linear triangle contains four mandatory
                # integers and one optional element index.

                if len(values) not in (4, 5):
                    raise ValueError(
                        f"Invalid triangle record {i}."
                    )

                element_type = int(values[0])

                # SU2/VTK type 5 corresponds to a linear
                # triangle with exactly three vertices.
                #
                # Quadrilaterals and other element types are
                # outside the scope of the current algorithm.

                if element_type != 5:
                    raise ValueError(
                        "Only linear triangular elements "
                        "(SU2 type 5) are supported."
                    )

                # Extract the three vertex indices.
                #
                # The element type and optional element
                # index are intentionally excluded.

                triangles[i] = [
                    int(values[1]),
                    int(values[2]),
                    int(values[3])
                ]

                # Validate the optional trailing element
                # index without storing it.
                #
                # Triangle numbering in our representation
                # is implicitly defined by the row index.

                if len(values) == 5:
                    int(values[4])

        # ====================================================
        # NPOIN — VERTEX COORDINATES
        # ====================================================

        # Example:
        #
        #     NPOIN= 4
        #     0.0 0.0 0
        #     1.0 0.0 1
        #     1.0 1.0 2
        #     0.0 1.0 3
        #
        # A two-dimensional point record contains:
        #
        #     x y [point_id]
        #
        # The global point numbering is implicitly defined
        # by the order of the records, starting at zero.
        #
        # The optional third value is a point index,
        # not a third spatial coordinate.

        elif key == "NPOIN":

            n_points = read_count(value)

            # Preallocate the vertex coordinate array.
            #
            # points[i, 0] = x coordinate of vertex i
            # points[i, 1] = y coordinate of vertex i

            points = np.empty(
                (n_points, 2),
                dtype=np.float64
            )

            for i in range(n_points):

                values = read_line().split()

                # Accept records with or without the optional
                # explicit point identifier.

                if len(values) not in (2, 3):
                    raise ValueError(
                        f"Invalid point record {i}."
                    )

                # Preserve the original coordinates and
                # the original vertex ordering.

                points[i] = [
                    float(values[0]),
                    float(values[1])
                ]

                # If an explicit point ID is provided,
                # verify that it agrees with the implicit
                # zero-based numbering.
                #
                # This implementation is intentionally limited
                # to complete, unpartitioned SU2 meshes.
                #
                # An inconsistent explicit ID may indicate
                # an unsupported numbering convention.

                if len(values) == 3:

                    point_id = int(values[2])

                    if point_id != i:
                        raise ValueError(
                            "Explicit point index does not match "
                            f"the implicit SU2 index: "
                            f"{point_id} != {i}."
                        )

        # ====================================================
        # NMARK — BOUNDARY MARKERS
        # ====================================================

        # Example:
        #
        #     NMARK= 2
        #
        #     MARKER_TAG= wall
        #     MARKER_ELEMS= 1
        #     3 0 1
        #
        #     MARKER_TAG= farfield
        #     MARKER_ELEMS= 3
        #     3 1 2
        #     3 2 3
        #     3 3 0
        #
        # NMARK specifies the number of boundary groups.
        #
        # Each boundary group contains:
        #
        #     MARKER_TAG   : physical boundary name
        #     MARKER_ELEMS : number of boundary elements
        #
        # A linear boundary edge has the format:
        #
        #     3 vertex_1 vertex_2
        #
        # Type 3 corresponds to a line segment.

        elif key == "NMARK":

            n_markers = read_count(value)

            for _ in range(n_markers):

                # Read the physical boundary name.

                _, tag = read_keyword(
                    read_line(),
                    "MARKER_TAG"
                )

                if not tag:
                    raise ValueError(
                        "Empty boundary marker."
                    )

                # Read the number of edges associated with
                # the current boundary marker.

                _, count = read_keyword(
                    read_line(),
                    "MARKER_ELEMS"
                )

                n_edges = read_count(count)

                for _ in range(n_edges):

                    values = read_line().split()

                    # Each linear boundary edge must contain
                    # its element type and two endpoint IDs.

                    if len(values) != 3:
                        raise ValueError(
                            "Invalid boundary edge record."
                        )

                    element_type = int(values[0])

                    if element_type != 3:
                        raise ValueError(
                            "Only linear boundary edges "
                            "(SU2 type 3) are supported."
                        )

                    # Store the connectivity of the edge.

                    boundary_edges.append([
                        int(values[1]),
                        int(values[2])
                    ])

                    # Store the associated marker name.
                    #
                    # Repeating the marker name for every
                    # edge creates a direct correspondence:
                    #
                    # boundary_tags[i] belongs to
                    # boundary_edges[i].
                    #
                    # This avoids storing nested dictionaries
                    # or custom data structures.

                    boundary_tags.append(tag)

        else:

            # Reject unknown sections instead of silently
            # skipping potentially important information.
            #
            # This includes unsupported multizone or
            # partitioned mesh conventions.

            raise ValueError(
                f"Unsupported SU2 mesh section: {key}"
            )

    # --------------------------------------------------------
    # 1.4. Convert boundary lists to NumPy arrays
    # --------------------------------------------------------

    # The reshape operation guarantees that the connectivity
    # array has two columns, even if no boundary edges exist.
    #
    # Expected shape:
    #
    #     (N_boundary_edges, 2)

    boundary_edges = np.asarray(
        boundary_edges,
        dtype=np.int64
    ).reshape(-1, 2)

    # Store the boundary names as NumPy Unicode strings.
    #
    # Avoid dtype=object because object arrays require pickle
    # when loading an NPZ file.

    boundary_tags = np.asarray(
        boundary_tags,
        dtype=str
    )

    # --------------------------------------------------------
    # 1.5. Validate the required mesh sections
    # --------------------------------------------------------

    required_sections = {
        "NDIME",
        "NELEM",
        "NPOIN",
        "NMARK"
    }

    missing_sections = required_sections - sections

    if missing_sections:
        raise ValueError(
            f"Missing mesh sections: {missing_sections}"
        )

    if len(points) == 0 or len(triangles) == 0:
        raise ValueError("The mesh is empty.")

    # --------------------------------------------------------
    # 1.6. Validate the vertex coordinates
    # --------------------------------------------------------

    # Reject NaN and infinite coordinates.
    #
    # These values would make subsequent geometric operations
    # such as gradient reconstruction undefined.

    if not np.all(np.isfinite(points)):
        raise ValueError(
            "Non-finite mesh coordinates."
        )

    n_points = len(points)

    # --------------------------------------------------------
    # 1.7. Validate the connectivity indices
    # --------------------------------------------------------

    # Every connectivity index must refer to an existing
    # vertex in the points array.
    #
    # Valid index range:
    #
    #     0 <= vertex_id < N_vertices
    #
    # Negative indices must also be rejected because NumPy
    # otherwise interprets them as indexing from the end
    # of an array.

    if (
        np.any(triangles < 0)
        or np.any(triangles >= n_points)
    ):
        raise ValueError(
            "Invalid triangle connectivity."
        )

    if (
        np.any(boundary_edges < 0)
        or np.any(boundary_edges >= n_points)
    ):
        raise ValueError(
            "Invalid boundary connectivity."
        )

    # --------------------------------------------------------
    # 1.8. Check for degenerate triangles
    # --------------------------------------------------------

    # Advanced NumPy indexing retrieves the coordinates of
    # the three vertices of every triangle.
    #
    # points.shape:
    #     (N_vertices, 2)
    #
    # triangles.shape:
    #     (N_triangles, 3)
    #
    # triangle_points.shape:
    #     (N_triangles, 3, 2)

    triangle_points = points[triangles]

    # Define two edges originating from the first vertex
    # of each triangle:
    #
    #     edge_1 = x_j - x_i
    #     edge_2 = x_l - x_i
    #
    # Both arrays have shape (N_triangles, 2).

    edge_1 = (
        triangle_points[:, 1]
        - triangle_points[:, 0]
    )

    edge_2 = (
        triangle_points[:, 2]
        - triangle_points[:, 0]
    )

    # Compute twice the signed area of every triangle
    # using the two-dimensional cross product:
    #
    # 2 * signed_area =
    #     edge_1_x * edge_2_y
    #     - edge_1_y * edge_2_x
    #
    # A positive result indicates counterclockwise ordering.
    # A negative result indicates clockwise ordering.
    # A zero result indicates a degenerate triangle.
    #
    # The signed area is used only for validation.
    # The original triangle orientation is preserved.

    double_area = (
        edge_1[:, 0] * edge_2[:, 1]
        - edge_1[:, 1] * edge_2[:, 0]
    )

    if not np.all(np.isfinite(double_area)):
        raise ValueError(
            "Non-finite triangle areas."
        )

    if np.any(double_area == 0.0):
        raise ValueError(
            "Degenerate triangles detected."
        )

    # --------------------------------------------------------
    # 1.9. Check boundary edge connectivity
    # --------------------------------------------------------

    # A boundary edge must connect two distinct vertices.
    #
    # An edge of the form [i, i] is invalid because it has
    # identical endpoint indices.

    if np.any(
        boundary_edges[:, 0] == boundary_edges[:, 1]
    ):
        raise ValueError(
            "Degenerate boundary edges detected."
        )

    # Return the original mesh represented exclusively
    # by numerical arrays.

    return (
        points,
        triangles,
        boundary_edges,
        boundary_tags
    )

def identify_boundary_nodes(boundary_edges, boundary_tags):
    """
    Identify protected vertices on the airfoil and farfield.

    Outlet vertices are excluded from the protected set,
    unless they also belong to a protected boundary.
    """

    # Select edges belonging to protected boundary markers.
    protected_markers = ("airfoil", "farfield")

    mask = np.isin(boundary_tags, protected_markers)

    # Extract unique vertex indices from protected edges.
    return np.unique(boundary_edges[mask].ravel())


# ============================================================
# 2. READ THE SU2 ASCII RESTART FILE
# ============================================================

def read_restart_field(restart_path, field_name):
    """
    Read a selected nodal scalar field from a SU2 restart CSV.

    Example CSV structure:

        "PointID","x","y","Density","Pressure","Mach"
        2,1.0,1.0,1.40,110000.0,0.60
        0,0.0,0.0,1.20,100000.0,0.80
        1,1.0,0.0,1.30,105000.0,0.95

    The CSV row ordering is not assumed to match the mesh
    vertex ordering.

    Parameters
    ----------
    restart_path : str or Path
        Path to the SU2 ASCII restart file.

    field_name : str
        Exact name of the scalar field to extract.

    Returns
    -------
    point_ids : ndarray of int64, shape (N_records,)
        Global SU2 point identifiers in CSV row order.

    coordinates : ndarray of float64, shape (N_records, 2)
        Restart coordinates in CSV row order.
        Returns None if coordinate columns are absent.

    values : ndarray of float64, shape (N_records,)
        Selected scalar values in CSV row order.

    Notes
    -----
    The scalar values are not reordered by this function.
    Mesh-to-solution matching is performed separately.
    """

    # Temporary lists are used because the number of nodal
    # records has not yet been established.

    point_ids = []
    coordinates = []
    values = []

    # Remove accidental whitespace from the requested name.

    field_name = field_name.strip()

    def clean_row(row):
        """
        Remove surrounding whitespace and an optional
        trailing CSV delimiter.

        Example:

            ["PointID", " x", " y", "Mach", ""]

        becomes:

            ["PointID", "x", "y", "Mach"]
        """

        row = [
            value.strip()
            for value in row
        ]

        if row and row[-1] == "":
            row.pop()

        return row

    # --------------------------------------------------------
    # 2.1. Open the restart file
    # --------------------------------------------------------

    # utf-8-sig also supports files beginning with a UTF-8
    # byte order mark.
    #
    # newline="" is recommended when using Python's csv module.

    with open(
        restart_path,
        "r",
        encoding="utf-8-sig",
        newline=""
    ) as file:

        # The csv module automatically handles quoted fields,
        # comma separation and standard CSV escaping.

        reader = csv.reader(
            file,
            skipinitialspace=True
        )

        # ----------------------------------------------------
        # 2.2. Read the CSV header
        # ----------------------------------------------------

        # The first record contains the column names.
        #
        # Example:
        #
        # PointID, x, y, Density, Pressure, Mach
        #
        # Column indices are determined dynamically, because
        # their order depends on the solver and output settings.

        header = next(reader, None)

        if header is None:
            raise ValueError(
                "The restart file is empty."
            )

        header = clean_row(header)

        # Duplicate column names would make field selection
        # ambiguous.

        if len(set(header)) != len(header):
            raise ValueError(
                "Duplicate CSV column names."
            )

        # ----------------------------------------------------
        # 2.3. Locate the PointID column
        # ----------------------------------------------------

        # PointID provides the direct association between
        # the restart solution and the original mesh.
        #
        # It must not be confused with the CSV row number.

        if "PointID" not in header:
            raise ValueError(
                "PointID column not found."
            )

        id_column = header.index("PointID")

        # ----------------------------------------------------
        # 2.4. Locate the requested scalar field
        # ----------------------------------------------------

        # The program is independent of the physical nature
        # of the scalar field.
        #
        # It can extract Mach, Density, Pressure, Temperature
        # or another available nodal quantity.
        #
        # No missing physical quantity is reconstructed.

        if field_name not in header:
            raise ValueError(
                f"Field '{field_name}' not found.\n"
                f"Available fields: {', '.join(header)}"
            )

        field_column = header.index(field_name)

        # ----------------------------------------------------
        # 2.5. Locate the coordinate columns
        # ----------------------------------------------------

        # Restart coordinates provide an independent check
        # that the solution belongs to the supplied mesh.
        #
        # They are not required for direct PointID matching,
        # but they are highly useful for detecting mismatches.

        has_x = "x" in header
        has_y = "y" in header

        # A two-dimensional coordinate record must provide
        # both x and y, or neither.

        if has_x != has_y:
            raise ValueError(
                "Incomplete coordinate columns."
            )

        has_coordinates = has_x and has_y

        if has_coordinates:

            x_column = header.index("x")
            y_column = header.index("y")

        # ----------------------------------------------------
        # 2.6. Read the nodal solution records
        # ----------------------------------------------------

        for row in reader:

            # Ignore completely empty CSV records.

            if (
                not row
                or not any(value.strip() for value in row)
            ):
                continue

            row = clean_row(row)

            # Every nodal record must contain the same
            # number of columns as the header.
            #
            # Unexpected records are rejected rather than
            # silently interpreted as solution data.

            if len(row) != len(header):
                raise ValueError(
                    f"Invalid CSV record at line "
                    f"{reader.line_num}."
                )

            # Extract the global point identifier.

            point_id = int(row[id_column])

            # Extract the selected physical scalar value.

            field_value = float(row[field_column])

            # Preserve the CSV row ordering for now.
            #
            # For example:
            #
            # point_ids = [2, 0, 1]
            # values    = [0.60, 0.80, 0.95]
            #
            # These arrays will be matched to the mesh in
            # the next processing stage.

            point_ids.append(point_id)
            values.append(field_value)

            # Extract the coordinates associated with the
            # same CSV record, if they are available.

            if has_coordinates:

                coordinates.append([
                    float(row[x_column]),
                    float(row[y_column])
                ])

    # --------------------------------------------------------
    # 2.7. Convert the extracted data to NumPy arrays
    # --------------------------------------------------------

    point_ids = np.asarray(
        point_ids,
        dtype=np.int64
    )

    values = np.asarray(
        values,
        dtype=np.float64
    )

    if has_coordinates:

        coordinates = np.asarray(
            coordinates,
            dtype=np.float64
        ).reshape(-1, 2)

    else:

        # None explicitly indicates that coordinate-based
        # verification is unavailable.

        coordinates = None

    # --------------------------------------------------------
    # 2.8. Validate the scalar field
    # --------------------------------------------------------

    # NaN and infinite values would propagate through the
    # future gradient and Hessian reconstruction.
    #
    # Reject them before constructing the adaptation input.

    if not np.all(np.isfinite(values)):
        raise ValueError(
            "Non-finite scalar field values."
        )

    return point_ids, coordinates, values

# ============================================================
# 3. MATCH THE MESH AND SCALAR FIELD
# ============================================================

def match_scalar_field(
    points,
    point_ids,
    coordinates,
    values
):
    """
    Associate each extracted scalar value with its mesh vertex.

    The restart file may contain the vertices in a different
    order from the mesh file.

    Example
    -------
    Mesh vertex ordering:

        points[0]
        points[1]
        points[2]

    Restart data:

        point_ids = [2, 0, 1]
        values    = [0.60, 0.80, 0.95]

    Correct output:

        scalar = [0.80, 0.95, 0.60]

    Parameters
    ----------
    points : ndarray, shape (N_vertices, 2)
        Original mesh coordinates.

    point_ids : ndarray, shape (N_records,)
        Global point identifiers in CSV row order.

    coordinates : ndarray or None
        Restart coordinates in CSV row order.

    values : ndarray, shape (N_records,)
        Selected scalar values in CSV row order.

    Returns
    -------
    scalar : ndarray of float64, shape (N_vertices,)
        Scalar field ordered consistently with the mesh.

    Fundamental invariant:

        scalar[i] = physical scalar value at points[i]
    """

    n_points = len(points)

    # --------------------------------------------------------
    # 3.1. Verify the number of nodal records
    # --------------------------------------------------------

    # A complete nodal restart must contain exactly one
    # solution record for every mesh vertex.

    if len(point_ids) != n_points:

        raise ValueError(
            "Mesh and solution point counts differ: "
            f"{n_points} != {len(point_ids)}."
        )

    # --------------------------------------------------------
    # 3.2. Verify the point identifier range
    # --------------------------------------------------------

    # Every restart PointID must reference an existing
    # mesh vertex.
    #
    # Valid IDs range from zero to N_vertices - 1.

    if (
        np.any(point_ids < 0)
        or np.any(point_ids >= n_points)
    ):
        raise ValueError(
            "PointID outside the valid mesh range."
        )

    # --------------------------------------------------------
    # 3.3. Detect missing or duplicate point identifiers
    # --------------------------------------------------------

    # np.bincount counts how many times each integer occurs.
    #
    # Example:
    #
    # point_ids   = [2, 0, 1]
    # occurrences = [1, 1, 1]
    #
    # This is valid because every point occurs once.
    #
    # Invalid example:
    #
    # point_ids   = [2, 0, 0]
    # occurrences = [2, 0, 1]
    #
    # Point 0 is duplicated and point 1 is missing.
    #
    # Checking only the number of CSV records would not
    # detect this error.

    occurrences = np.bincount(
        point_ids,
        minlength=n_points
    )

    if np.any(occurrences != 1):
        raise ValueError(
            "Missing or duplicate PointIDs."
        )

    # --------------------------------------------------------
    # 3.4. Verify the restart coordinates
    # --------------------------------------------------------

    # When coordinates are available in the restart file,
    # verify that every PointID references the correct
    # physical location.
    #
    # This detects situations where two different meshes
    # have the same number of vertices and point IDs.

    if coordinates is not None:

        if not np.all(np.isfinite(coordinates)):
            raise ValueError(
                "Non-finite restart coordinates."
            )

        # ----------------------------------------------------
        # Determine the characteristic domain size
        # ----------------------------------------------------

        # np.ptp calculates the coordinate range:
        #
        #     max(coordinate) - min(coordinate)
        #
        # Applied along axis 0, it returns:
        #
        #     [x_max - x_min, y_max - y_min]
        #
        # Taking the maximum provides a characteristic
        # length scale for the computational domain.

        domain_size = np.max(
            np.ptp(points, axis=0)
        )

        # Determine the maximum absolute coordinate value.
        #
        # This is useful when estimating floating-point
        # round-off for domains located far from the origin.

        coordinate_scale = np.max(
            np.abs(points)
        )

        # ----------------------------------------------------
        # Define a geometric comparison tolerance
        # ----------------------------------------------------

        # Exact floating-point equality is not appropriate
        # because the mesh and restart files may use
        # different output precision.
        #
        # The tolerance contains:
        #
        # 1. A domain-relative geometric tolerance.
        # 2. A small floating-point round-off allowance.
        #
        # The factor 1e-9 is an implementation choice and
        # should be checked against actual SU2 output precision.

        tolerance = (
            1e-9 * domain_size
            + 16
            * np.finfo(np.float64).eps
            * coordinate_scale
        )

        # ----------------------------------------------------
        # Match coordinates using the global PointIDs
        # ----------------------------------------------------

        # Advanced indexing retrieves the mesh coordinates
        # in the same order as the restart records.
        #
        # Example:
        #
        # point_ids = [2, 0, 1]
        #
        # points[point_ids] =
        #
        #     [
        #         points[2],
        #         points[0],
        #         points[1]
        #     ]
        #
        # Both coordinate arrays can now be compared
        # directly, component by component.

        mesh_coordinates = points[point_ids]

        coordinate_errors = np.abs(
            coordinates - mesh_coordinates
        )

        # coordinate_errors has shape (N_vertices, 2).
        #
        # Its maximum identifies the largest absolute
        # coordinate discrepancy.

        max_error = np.max(coordinate_errors)

        # ----------------------------------------------------
        # Report an inconsistent mesh-to-solution association
        # ----------------------------------------------------

        if max_error > tolerance:

            # Find the CSV record and coordinate component
            # responsible for the maximum discrepancy.

            row, component = np.unravel_index(
                np.argmax(coordinate_errors),
                coordinate_errors.shape
            )

            coordinate_name = ("x", "y")[component]

            raise ValueError(
                f"Coordinate mismatch at "
                f"PointID {point_ids[row]}.\n"
                f"Coordinate component: {coordinate_name}\n"
                f"Mesh coordinates: "
                f"{points[point_ids[row]]}\n"
                f"Restart coordinates: "
                f"{coordinates[row]}\n"
                f"Maximum error: {max_error:.6e}\n"
                f"Tolerance: {tolerance:.6e}"
            )

    else:

        # PointID matching remains possible without restart
        # coordinates.
        #
        # However, the geometric identity of the two files
        # cannot be independently verified.

        print(
            "Warning: no restart coordinates available. "
            "Matching uses PointID only."
        )

    # --------------------------------------------------------
    # 3.5. Reorder the scalar field
    # --------------------------------------------------------

    # Allocate one scalar value per mesh vertex.
    #
    # np.empty does not initialise its entries.
    #
    # This is safe because the previous checks guarantee
    # that every valid point ID occurs exactly once.

    scalar = np.empty(
        n_points,
        dtype=np.float64
    )

    # This is the fundamental mesh-to-solution matching
    # operation.
    #
    # Example:
    #
    # point_ids = [2, 0, 1]
    # values    = [0.60, 0.80, 0.95]
    #
    # The vectorised assignment performs:
    #
    # scalar[2] = 0.60
    # scalar[0] = 0.80
    # scalar[1] = 0.95
    #
    # Final result:
    #
    # scalar = [0.80, 0.95, 0.60]
    #
    # No sorting, spatial search or interpolation is needed.

    scalar[point_ids] = values

    # The scalar field and the mesh coordinates now use
    # exactly the same vertex indexing convention.

    return scalar

# ============================================================
# 4. EXTRACT AND EXPORT THE DATA
# ============================================================

def extract_solution(
    mesh_path,
    restart_path,
    field_name,
    output_path
):
    """
    Convert a SU2 mesh and an associated nodal scalar field
    into a single NumPy NPZ file.

    This function coordinates the complete extraction process:

        1. Read the SU2 mesh.
        2. Read the requested restart field.
        3. Match the scalar values to the mesh vertices.
        4. Export the validated arrays.

    Parameters
    ----------
    mesh_path : str or Path
        Input SU2 ASCII mesh file.

    restart_path : str or Path
        Input SU2 ASCII restart CSV file.

    field_name : str
        Name of the nodal scalar field to extract.

    output_path : str or Path
        Destination NPZ file.

    Example
    -------
    extract_solution(
        "mesh.su2",
        "restart_flow.csv",
        "Mach",
        "adaptation_input.npz"
    )
    """

    # Convert the input paths to pathlib Path objects.
    #
    # This allows consistent path manipulation across
    # macOS, Linux and Windows.

    mesh_path = Path(mesh_path)
    restart_path = Path(restart_path)
    output_path = Path(output_path)

    field_name = field_name.strip()

    # --------------------------------------------------------
    # 4.1. Validate the output path
    # --------------------------------------------------------

    if output_path.suffix.lower() != ".npz":
        raise ValueError(
            "The output file must have a .npz extension."
        )

    # resolve() converts paths into absolute paths and
    # resolves symbolic links.
    #
    # This prevents the output path from accidentally
    # referring to one of the original input files.

    if output_path.resolve() in (
        mesh_path.resolve(),
        restart_path.resolve()
    ):
        raise ValueError(
            "The output would overwrite an input file."
        )

    # --------------------------------------------------------
    # 4.2. Read the original mesh
    # --------------------------------------------------------

    (
        points,
        triangles,
        boundary_edges,
        boundary_tags
    ) = read_su2_mesh(mesh_path)

    # Collect the vertex IDs directly from all declared SU2 markers.
    # These node IDs will be protected during mesh adaptation.

    boundary_nodes = identify_boundary_nodes(boundary_edges, boundary_tags)

    # At this stage, the following arrays are available:
    #
    # points         : (N_vertices, 2)
    # triangles      : (N_triangles, 3)
    # boundary_edges : (N_boundary_edges, 2)
    # boundary_tags  : (N_boundary_edges,)
    # boundary_nodes : (N_boundary_nodes,)

    # --------------------------------------------------------
    # 4.3. Extract the requested physical scalar field
    # --------------------------------------------------------

    point_ids, coordinates, values = read_restart_field(
        restart_path,
        field_name
    )

    # The extracted scalar values still follow the original
    # restart CSV row ordering.

    # --------------------------------------------------------
    # 4.4. Associate the scalar field with the mesh
    # --------------------------------------------------------

    scalar = match_scalar_field(
        points,
        point_ids,
        coordinates,
        values
    )

    # The fundamental correspondence is now established:
    #
    # points[i] = coordinates of mesh vertex i
    # scalar[i] = physical scalar value at mesh vertex i

    # --------------------------------------------------------
    # 4.5. Export the NPZ file
    # --------------------------------------------------------

    # An NPZ file is a ZIP archive containing individual
    # NumPy arrays stored in the NPY format.
    #
    # Unlike a CSV file, it preserves the shape and numerical
    # data type of every array.
    #
    # The compressed version reduces disk storage without
    # introducing any numerical approximation.
    #
    # Exactly seven arrays are stored.
    #
    # No custom objects, dictionaries or pickle-based
    # representations are required.

    with open(output_path, "wb") as file:

        np.savez_compressed(
            file,

            # Cartesian vertex coordinates.
            points=points,

            # Triangular element connectivity.
            triangles=triangles,

            # Nodal scalar field aligned with the vertices.
            scalar=scalar,

            # Boundary edge connectivity.
            boundary_edges=boundary_edges,

            # Physical marker associated with every edge.
            boundary_tags=boundary_tags,

            # Fixed vertex indices on all physical boundaries.
            boundary_nodes=boundary_nodes,

            # Zero-dimensional Unicode array containing
            # the name of the extracted scalar field.
            field_name=np.array(
                field_name,
                dtype=str
            )
        )

    # --------------------------------------------------------
    # 4.6. Print the extraction summary
    # --------------------------------------------------------

    # These statistics provide a quick overview of the
    # generated adaptation input.
    #
    # They are not mesh adaptation calculations.

    print("\nSU2 solution extracted successfully.\n")

    print(f"Scalar field:     {field_name}")
    print(f"Vertices:         {len(points)}")
    print(f"Triangles:        {len(triangles)}")
    print(f"Boundary edges:   {len(boundary_edges)}")
    print(f"Boundary nodes:   {len(boundary_nodes)}")

    # np.unique returns the distinct marker names.

    print(
        f"Boundary markers: "
        f"{', '.join(np.unique(boundary_tags))}"
    )

    print(
        f"Scalar minimum:   {np.min(scalar):.8e}"
    )

    print(
        f"Scalar maximum:   {np.max(scalar):.8e}"
    )

    print(f"Output:           {output_path}")

# ============================================================
# 5. COMMAND-LINE INTERFACE
# ============================================================

def main():
    """
    Provide a command-line interface to extract_solution().

    Usage
    -----
    python extract_solution.py \\
        mesh.su2 \\
        restart_flow.csv \\
        Mach \\
        adaptation_input.npz
    """

    # argparse provides standard command-line argument
    # handling and automatically generates a help message.
    #
    # All four input arguments are positional and mandatory.

    parser = argparse.ArgumentParser(
        description=(
            "Extract a SU2 triangular mesh and a nodal "
            "scalar field into a NumPy NPZ file."
        )
    )

    parser.add_argument(
        "mesh_path",
        help="Path to the input SU2 mesh."
    )

    parser.add_argument(
        "restart_path",
        help="Path to the SU2 ASCII restart CSV."
    )

    parser.add_argument(
        "field_name",
        help="Name of the scalar field to extract."
    )

    parser.add_argument(
        "output_path",
        help="Path to the output NPZ file."
    )

    # Parse the arguments provided by the user.

    args = parser.parse_args()

    # Execute the extraction.
    #
    # Expected input/output and parsing errors are converted
    # into concise command-line error messages.
    #
    # Exit status:
    #
    # 0 : successful execution
    # 1 : extraction failure

    try:

        extract_solution(
            args.mesh_path,
            args.restart_path,
            args.field_name,
            args.output_path
        )

    except (
        OSError,
        ValueError,
        OverflowError,
        csv.Error
    ) as error:

        parser.exit(
            1,
            f"Error: {error}\n"
        )

# ============================================================
# 6. PROGRAM ENTRY POINT
# ============================================================

# This condition distinguishes two situations:
#
# 1. Direct execution:
#
#       python extract_solution.py ...
#
#    The command-line interface is executed.
#
# 2. Import from another Python script:
#
#       from extract_solution import extract_solution
#
#    Only the function definitions are loaded.
#    No file is read or written automatically.
#
# This allows the extraction function to be reused directly
# inside the future automated CFD adaptation workflow.

if __name__ == "__main__":
    main()
