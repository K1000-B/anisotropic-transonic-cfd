"""Orchestrate steady, single-zone SU2 solve/extract/adapt steps on the host.

Implementation-only mission: this source has not been executed or tested.
The elementary programs remain authoritative for extraction, adaptation and
Docker execution. Importing this module never starts a solver or prompts.
"""

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXPECTED_ENVIRONMENT = PROJECT_ROOT / ".PROJECT_CFDvenv"

# Check the selected host interpreter before importing scientific dependencies.
# Do not resolve sys.executable: the venv executable can be a symlink.
if (Path(sys.prefix).absolute() != EXPECTED_ENVIRONMENT
        or Path(sys.executable).absolute().parent != EXPECTED_ENVIRONMENT / "bin"):
    raise RuntimeError("Use this checkout's .PROJECT_CFDvenv/bin/python.")

import numpy as np

# Numbered stage directories are intentionally kept in their existing layout.
sys.path.insert(0, str(PROJECT_ROOT / "02_mesh_adaptation"))
from extract_solution import (
    extract_solution, match_scalar_field, read_restart_field, read_su2_mesh,
)
from mesh_adaptation import bounded_metric, mesh_adaptation, triangle_geometry
from run_simulation import (
    CONTAINER_NAME, host_to_container_path, inspect_container, run_simulation,
)


CYCLE_KEYS = ("SUBITER", "BUDGET", "FLOW_ITER", "HMINS", "P_NORMS", "ARMAXS")
GLOBAL_KEYS = ("SENSOR", "CHORD", "HMAX", "EPSILON", "COLLAPSE_THRESHOLD",
               "SPLIT_THRESHOLD", "REMESH_MAX_ITER")
HISTORY_FIELDS = (
    "solver case_name run_id global_step cycle subiteration sensor target_nodes "
    "actual_nodes actual_triangles target_complexity h_min h_max p epsilon "
    "max_aspect_ratio collapse_threshold split_threshold remesh_max_iterations "
    "flow_iteration_limit remeshing_status extraction_runtime_seconds "
    "adaptation_runtime_seconds plot_runtime_seconds cfd_runtime_seconds "
    "total_runtime_seconds mesh_path restart_path vtu_path "
    "bounded_metric_complexity complexity_relative_error metric_quality_min "
    "metric_quality_median physical_quality_min physical_quality_median "
    "metric_edge_rms short_edge_count long_edge_count frozen_boundary_violations "
    "split_count collapse_count swap_count relocation_count "
    "actual_flow_iterations cfd_convergence_status final_history_values"
).split()

# Filenames that SU2 writes, rather than reads. All are confined to a solve's
# fresh directory; their original contents/formats remain otherwise unchanged.
OUTPUT_PATH_KEYS = {
    "CONV_FILENAME", "RESTART_FILENAME", "RESTART_ADJ_FILENAME", "VOLUME_FILENAME",
    "VOLUME_ADJ_FILENAME", "SURFACE_FILENAME", "SURFACE_ADJ_FILENAME",
    "SURFACE_SENS_FILENAME", "VOLUME_SENS_FILENAME", "MESH_OUT_FILENAME",
    "BREAKDOWN_FILENAME", "VALUE_OBJFUNC_FILENAME", "GRAD_OBJFUNC_FILENAME",
    "HESS_OBJFUNC_FILENAME",
}
UNUSED_RESTART_KEYS = {"SOLUTION_FILENAME", "SOLUTION_ADJ_FILENAME"}


# ---------------------------------------------------------------------------
# Input parsing: no shell evaluation and no eval().
# ---------------------------------------------------------------------------
def resolve_input_file(value, allowed_extensions):
    """Accept normal, quoted and terminal-escaped host paths without a shell."""
    if isinstance(value, Path):
        path = value.expanduser()
    else:
        raw = str(value).strip()
        # Prefer an existing literal path, including a literal backslash/name.
        path = Path(raw).expanduser()
        if not path.is_file():
            tokens = shlex.split(raw)
            if len(tokens) != 1:
                raise ValueError("Supply exactly one file path.")
            path = Path(tokens[0]).expanduser()
    path = path.resolve(strict=True)
    if not path.is_file() or path.suffix.lower() not in allowed_extensions:
        raise ValueError(f"Expected a regular file with extension {allowed_extensions}: {path}")
    with path.open("rb") as stream:
        stream.read(1)  # Also check actual read access.
    return path


def request_input_file(prompt, allowed_extensions, validator):
    """Retry one input until both its path and its contents are acceptable."""
    while True:
        try:
            path = resolve_input_file(input(prompt + "\n> "), allowed_extensions)
            value = validator(path)
            return path, value
        except (OSError, ValueError, OverflowError, csv.Error) as error:
            print(f"Invalid input: {error}", flush=True)


def read_assignments(config_path):
    """Read single-line SU2-style assignments, retaining the original text."""
    text = Path(config_path).read_text(encoding="utf-8-sig")
    values = {}
    for number, original in enumerate(text.splitlines(), 1):
        line = original.split("%", 1)[0].strip()
        if not line:
            continue
        key, separator, value = line.partition("=")
        key, value = key.strip(), value.strip()
        if not separator or not re.fullmatch(r"[A-Z][A-Z0-9_]*", key) or not value:
            raise ValueError(f"Malformed assignment at {config_path}:{number}.")
        if key in values:
            raise ValueError(f"Duplicate active parameter {key} at line {number}.")
        values[key] = value
    return text, values


def value_tokens(value):
    """Split one scalar/list value, respecting quoted strings and commas."""
    if value.startswith("("):
        if not value.endswith(")"):
            raise ValueError(f"Unclosed parameter list: {value}")
        value = value[1:-1].strip()
    elif "(" in value or ")" in value:
        raise ValueError(f"Malformed parameter value: {value}")
    tokens = next(csv.reader([value], skipinitialspace=True, strict=True))
    tokens = [token.strip() for token in tokens]
    if not tokens or any(not token or "(" in token or ")" in token for token in tokens):
        raise ValueError(f"Empty or malformed parameter value: {value}")
    return tokens


def read_adaptation_config(config_path):
    """Parse the complete, case-owned adaptation controls into Python values."""
    _, raw = read_assignments(config_path)
    expected = set(CYCLE_KEYS + GLOBAL_KEYS)
    missing, unknown = expected - raw.keys(), raw.keys() - expected
    if missing or unknown:
        raise ValueError(f"Missing parameters: {sorted(missing)}; unknown: {sorted(unknown)}.")
    parameters = {}
    for key, value in raw.items():
        tokens = value_tokens(value)
        if key == "SENSOR":
            if len(tokens) != 1:
                raise ValueError("SENSOR must name exactly one scalar field.")
            parameters[key] = tokens[0]
            continue
        numbers = []
        for token in tokens:
            # Integer syntax stays integer; 2.0 is converted only where allowed.
            numbers.append(int(token) if re.fullmatch(r"[+-]?\d+", token) else float(token))
        if key in CYCLE_KEYS:
            if not value.startswith("("):
                raise ValueError(f"{key} must be a parenthesised cycle list.")
            parameters[key] = numbers
        else:
            if len(numbers) != 1 or value.startswith("("):
                raise ValueError(f"{key} must be a scalar.")
            parameters[key] = numbers[0]
    validate_adaptation_config(parameters)
    return parameters


def validate_adaptation_config(parameters):
    """Validate every numerical control without inventing default parameters."""
    count = len(parameters["SUBITER"])
    if count == 0 or any(len(parameters[key]) != count for key in CYCLE_KEYS):
        raise ValueError("Cycle lists must have the same nonzero length.")
    for key in CYCLE_KEYS + GLOBAL_KEYS:
        if key == "SENSOR":
            if not parameters[key].strip():
                raise ValueError("SENSOR cannot be empty.")
            continue
        values = parameters[key] if key in CYCLE_KEYS else [parameters[key]]
        if any(not math.isfinite(value) for value in values):
            raise ValueError(f"{key} contains a non-finite value.")
    for key in ("SUBITER", "BUDGET", "FLOW_ITER", "REMESH_MAX_ITER"):
        values = parameters[key] if key in CYCLE_KEYS else [parameters[key]]
        if any(type(value) is not int or value <= 0 for value in values):
            raise ValueError(f"{key} requires positive integers.")
    if any(value <= 0 or value != int(value) for value in parameters["P_NORMS"]):
        raise ValueError("P_NORMS requires positive integer-valued numbers.")
    parameters["P_NORMS"] = [int(value) for value in parameters["P_NORMS"]]
    if any(parameters[key] <= 0 for key in ("CHORD", "HMAX", "EPSILON")):
        raise ValueError("CHORD, HMAX and EPSILON must be positive.")
    if any(not 0 < value <= parameters["HMAX"] for value in parameters["HMINS"]):
        raise ValueError("Each HMIN must satisfy 0 < HMIN <= HMAX.")
    if any(value < 1 for value in parameters["ARMAXS"]):
        raise ValueError("Every ARMAX must be at least one.")
    if not 0 < parameters["COLLAPSE_THRESHOLD"] < 1 < parameters["SPLIT_THRESHOLD"]:
        raise ValueError("Require 0 < COLLAPSE_THRESHOLD < 1 < SPLIT_THRESHOLD.")
    lengths = [parameters["HMAX"] * parameters["CHORD"]] + [
        value * parameters["CHORD"] for value in parameters["HMINS"]
    ]
    if any(not math.isfinite(value) or value <= 0 for value in lengths):
        raise ValueError("Chord-scaled metric sizes must be finite and positive.")


def validate_mesh(mesh_path):
    """Reuse the existing triangular SU2 reader and geometry safeguards."""
    data = read_su2_mesh(mesh_path)
    _, areas = triangle_geometry(data[0], data[1].copy())
    if not np.isfinite(areas.sum()) or areas.sum() <= 0:
        raise ValueError("Invalid physical mesh area.")
    if not np.any(data[3] == "airfoil"):
        raise ValueError("An airfoil boundary marker is required for fixed plot limits.")
    return data


def detect_cfd_solver(config_path):
    """Accept only steady, stationary, single-zone direct EULER/RANS cases."""
    _, values = read_assignments(config_path)
    solver = values.get("SOLVER", "EULER").upper()
    if solver not in ("EULER", "RANS"):
        raise ValueError(f"Unsupported SOLVER: {solver}.")
    for key, default, accepted in (
        ("TIME_DOMAIN", "NO", {"NO"}),
        ("TIME_MARCHING", "NO", {"NO"}),
        ("MATH_PROBLEM", "DIRECT", {"DIRECT"}),
        ("MESH_FORMAT", "SU2", {"SU2"}),
        ("GRID_MOVEMENT", "NO", {"NO", "NONE"}),
        ("DEFORM_MESH", "NO", {"NO"}),
    ):
        if values.get(key, default).upper() not in accepted:
            raise ValueError(f"{key} is incompatible with this steady fixed-mesh workflow.")
    if "CONFIG_LIST" in values and values["CONFIG_LIST"].upper() != "NONE":
        raise ValueError("Multizone CONFIG_LIST is unsupported.")
    return solver


def check_metric_feasibility(mesh_data, parameters):
    """Reuse bounded_metric's own feasibility check before paying for CFD."""
    _, areas = triangle_geometry(mesh_data[0], mesh_data[1].copy())
    # Only total area and size bounds enter feasibility. A one-node identity
    # probe invokes the existing implementation; no sensor/metric is rebuilt.
    for cycle, budget in enumerate(parameters["BUDGET"]):
        h_min = parameters["HMINS"][cycle] * parameters["CHORD"]
        h_max = parameters["HMAX"] * parameters["CHORD"]
        bounded_metric(
            np.eye(2)[None, :, :], np.array([areas.sum()]),
            target_complexity(budget), h_min, h_max, parameters["ARMAXS"][cycle],
        )


def target_complexity(target_nodes):
    """Report the existing algorithm's indicative-budget conversion."""
    return float(np.sqrt(3.0) / 2.0 * target_nodes)


# ---------------------------------------------------------------------------
# Preserve inputs and resolve paths before any solver is launched.
# ---------------------------------------------------------------------------
def resolve_dependencies(config_path):
    """Resolve explicit file inputs relative to the original configuration.

    Scalar file-valued options are supported. Nested configurations, lists of
    files and unrecognised FILE options fail explicitly instead of losing paths.
    Dependencies that themselves refer to further files need a later review.
    """
    _, values = read_assignments(config_path)
    dependencies = {}
    for key, value in values.items():
        if key in OUTPUT_PATH_KEYS | UNUSED_RESTART_KEYS | {"MESH_FILENAME", "OUTPUT_FILES"}:
            continue
        if not ("FILENAME" in key or key.endswith(("_FILE", "_FILES"))):
            continue
        if value.upper() in ("NONE", "NO"):
            continue
        if key in {"CONFIG_LIST", "RUNTIME_FILENAME"} or key.endswith("_FILES"):
            raise ValueError(f"Unsupported external file mechanism: {key}.")
        tokens = value_tokens(value)
        if len(tokens) != 1 or value.startswith("("):
            raise ValueError(f"Only scalar input filenames are supported: {key}.")
        source = Path(tokens[0]).expanduser()
        if not source.is_absolute():
            source = config_path.parent / source
        source = source.resolve(strict=True)
        if not source.is_file():
            raise ValueError(f"Required input is not a regular file: {key}: {source}")
        with source.open("rb") as stream:
            stream.read(1)
        if source.suffix.lower() == ".cfg":
            raise ValueError(f"Nested configuration dependencies require review: {key}.")
        dependencies[key] = source
    return dependencies


def snapshot_dependencies(dependencies, directory):
    """Copy resolved files under safe relative names, retaining exact bytes."""
    rewritten = {}
    if dependencies:
        (directory / "inputs").mkdir()
    for number, (key, source) in enumerate(sorted(dependencies.items()), 1):
        suffix = source.suffix if re.fullmatch(r"\.[A-Za-z0-9]+", source.suffix) else ".dat"
        relative = Path("inputs") / f"input_{number:03d}{suffix}"
        shutil.copy2(source, directory / relative)
        rewritten[key] = relative.as_posix()
    return rewritten


def resolve_case_directory(project_root, solver, mesh_path, config_path,
                           adaptation_config_path, dependencies):
    """Reuse an existing category case or create a collision-free case copy."""
    category = project_root / ("04_euler_simulations" if solver == "EULER"
                               else "05_rans_simulations")
    relative = None
    try:
        relative = config_path.relative_to(category.resolve())
    except ValueError:
        pass
    if relative is not None and len(relative.parts) >= 2:
        if "results" in relative.parts[:-1]:
            raise ValueError("Select case source inputs, not a previous results directory.")
        return category / relative.parts[0]

    # Never merge external inputs into an unrelated pre-existing case.
    name = config_path.stem if config_path.parent == category.resolve() else config_path.parent.name
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("._-") or config_path.stem
    index = 0
    while True:
        case = category / (name if index == 0 else f"{name}_{index:03d}")
        try:
            case.mkdir()
            break
        except FileExistsError:
            index += 1
    shutil.copy2(mesh_path, case / "baseline.su2")
    shutil.copy2(adaptation_config_path, case / "adaptation.cfg")
    shutil.copy2(config_path, case / "reference_config.cfg")
    paths = snapshot_dependencies(dependencies, case)
    paths["MESH_FILENAME"] = "baseline.su2"
    (case / "config.cfg").write_text(rewrite_configuration(config_path, paths), encoding="utf-8")
    return case


def create_run_directory(case_directory):
    """Reserve the next run number atomically, including partially filled runs."""
    results = case_directory / "results"
    results.mkdir(exist_ok=True)
    numbers = [int(match.group(1)) for path in results.iterdir()
               if (match := re.fullmatch(r"run_(\d+)", path.name))]
    number = max(numbers, default=0) + 1
    while True:
        run = results / f"run_{number:03d}"
        try:
            run.mkdir()
            return run
        except FileExistsError:
            number += 1


def rewrite_configuration(reference_path, replacements):
    """Preserve physical settings/comments; replace only selected assignments."""
    text, _ = read_assignments(reference_path)
    remaining = dict(replacements)
    lines = []
    for line in text.splitlines():
        active, marker, comment = line.partition("%")
        key = active.partition("=")[0].strip()
        if "=" in active and key in remaining:
            line = f"{key}= {remaining.pop(key)}"
            if marker:
                line += "  %" + comment
        lines.append(line)
    lines.extend(f"{key}= {value}" for key, value in remaining.items())
    return "\n".join(lines) + "\n"


def prepare_cfd_config(reference_config_path, mesh_path, flow_iterations,
                       output_directory, solver, sensor, dependencies):
    """Write a retained local config, changing execution/output controls only."""
    _, values = read_assignments(reference_config_path)
    replacements = snapshot_dependencies(dependencies, output_directory)
    replacements.update({
        "SOLVER": solver, "MESH_FILENAME": mesh_path.name, "ITER": str(flow_iterations),
        "RESTART_SOL": "NO", "RESTART_ITER": "0", "READ_BINARY_RESTART": "NO",
        "WRT_RESTART_COMPACT": "NO", "TABULAR_FORMAT": "CSV",
        "RESTART_FILENAME": "restart_flow", "VOLUME_FILENAME": "flow",
        "CONV_FILENAME": "history", "WRT_RESTART_OVERWRITE": "YES",
        "WRT_VOLUME_OVERWRITE": "YES", "HISTORY_WRT_FREQ_INNER": "1",
    })
    # Preserve additional requested output formats and their frequency pairing.
    outputs = value_tokens(values.get("OUTPUT_FILES", "(RESTART, PARAVIEW, SURFACE_PARAVIEW)"))
    outputs = [item.upper() for item in outputs]
    if len(set(outputs)) != len(outputs) or "NONE" in outputs:
        raise ValueError("Duplicate/NONE OUTPUT_FILES entries require correction.")
    frequencies = None
    if "OUTPUT_WRT_FREQ" in values:
        frequencies = value_tokens(values["OUTPUT_WRT_FREQ"])
        if len(frequencies) != len(outputs) or any(not item.isdigit() or int(item) <= 0
                                                  for item in frequencies):
            raise ValueError("OUTPUT_WRT_FREQ must pair positive integers with OUTPUT_FILES.")
    for required in ("RESTART_ASCII", "PARAVIEW"):
        if required not in outputs:
            outputs.append(required)
            if frequencies is not None:
                frequencies.append(str(flow_iterations))
    replacements["OUTPUT_FILES"] = "(" + ", ".join(outputs) + ")"
    if frequencies is not None:
        replacements["OUTPUT_WRT_FREQ"] = "(" + ", ".join(frequencies) + ")"
    volume = value_tokens(values.get("VOLUME_OUTPUT", "(COORDINATES, SOLUTION, PRIMITIVE)"))
    # COMPACT deliberately suppresses required fields in some outputs; remove
    # that output-size directive while keeping every physical output group.
    volume = [item for item in volume if item.upper() not in {"COMPACT", "NONE"}]
    for group in ("COORDINATES", "SOLUTION", "PRIMITIVE"):
        if group not in [item.upper() for item in volume]:
            volume.append(group)
    # Generic sensors use the user's preserved output groups. CSV display
    # names are not necessarily SU2 identifiers: never guess that translation.
    if sensor.upper() not in {"MACH", "DENSITY", "PRESSURE", "TEMPERATURE"}:
        print(f"Sensor {sensor!r} requires a field supplied by the preserved VOLUME_OUTPUT groups.", flush=True)
    replacements["VOLUME_OUTPUT"] = "(" + ", ".join(volume) + ")"
    history = value_tokens(values.get("HISTORY_OUTPUT", "(ITER, RMS_RES)"))
    history = [item for item in history if item.upper() != "NONE"]
    if "ITER" not in [item.upper() for item in history]:
        history.append("ITER")
    replacements["HISTORY_OUTPUT"] = "(" + ", ".join(history) + ")"
    for key in OUTPUT_PATH_KEYS - replacements.keys():
        if key in values:
            basename = Path(value_tokens(values[key])[0]).name
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", basename) or basename in (".", ".."):
                raise ValueError(f"Unsafe output basename for {key}.")
            # Distinct auxiliary stems avoid collisions with mesh/config/flow
            # files and with other requested output formats.
            replacements[key] = "aux_" + key.lower() + Path(basename).suffix
    path = output_directory / "config.cfg"
    path.write_text(rewrite_configuration(reference_config_path, replacements), encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Output integrity, fixed plotting limits, progressive PVD and CSV recording.
# ---------------------------------------------------------------------------
def resolve_sensor_name(restart_path, sensor):
    """Find one case-insensitive match, then pass its exact spelling onward."""
    with restart_path.open(encoding="utf-8-sig", newline="") as stream:
        header = next(csv.reader(stream, skipinitialspace=True), [])
    header = [name.strip() for name in header]
    while header and not header[-1]:
        header.pop()
    matches = [name for name in header if name.casefold() == sensor.casefold()]
    if len(matches) != 1:
        raise ValueError(f"Sensor {sensor!r} has {len(matches)} matches. Available fields: {header}")
    return matches[0]


def require_output(path):
    """Reject absent/empty/unreadable outputs; a fresh solve has no stale files."""
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Required nonempty output missing: {path}")
    with path.open("rb") as stream:
        stream.read(1)


def inspect_history(history_path, log_path, flow_iterations):
    """Record final history values and explicit solver convergence messages.

    A limit or an early exit does not itself prove convergence. Keep an unknown
    status if the version-specific success message cannot be identified.
    """
    last = None
    with history_path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, skipinitialspace=True)
        header = [value.strip() for value in next(reader, [])]
        if not header or len(set(header)) != len(header):
            raise ValueError("Missing or duplicate SU2 history header.")
        for row in reader:
            if not row or not any(value.strip() for value in row):
                continue
            if len(row) != len(header):
                raise ValueError("Malformed SU2 convergence history row.")
            numbers = [float(value) for value in row]
            if not all(math.isfinite(value) for value in numbers):
                raise ValueError("Non-finite CFD history: refusing further adaptation.")
            last = dict(zip(header, numbers))
    if last is None:
        raise ValueError("SU2 history has no iteration data.")
    # CSinglezoneDriver.cpp v8.5.0 reports these two outcomes explicitly.
    converged, limited = False, False
    with log_path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, 1):
            if "Interrupt signal (" in line:
                raise ValueError("SU2 reported an interruption; this solve is incomplete.")
            converged |= "All convergence criteria satisfied." in line
            limited |= "Maximum number of iterations reached (ITER =" in line
            # SU2 describes freestream initialisation as "values at infinity".
            # Require a failure diagnostic or a delimited numerical value;
            # searching for bare inf/infinity in arbitrary prose is unsafe.
            failure = re.search(
                r"\bdiverged\b|\bdivergence detected\b|\bnon-finite\b|"
                r"\b(?:nan|inf(?:inity)?)\s+(?:detected|encountered)\b",
                line, re.IGNORECASE,
            )
            nonfinite_value = re.search(
                r"(?:^|[|=:])\s*[+-]?(?:nan|inf(?:inity)?)\s*(?=$|[|,;])",
                line, re.IGNORECASE,
            )
            if failure or nonfinite_value:
                raise ValueError(
                    "Solver log reports divergence/non-finite values at "
                    f"{log_path.name}:{line_number}: {line.strip()}"
                )
    iteration = last.get("Inner_Iter")
    if iteration is not None and (iteration < 0 or iteration != int(iteration)):
        raise ValueError("Invalid Inner_Iter history value.")
    actual = int(iteration) + 1 if iteration is not None else ""
    status = "converged_by_solver_report" if converged else (
        "iteration_limit_without_confirmed_convergence" if limited or (
            actual != "" and actual >= flow_iterations) else "unconfirmed"
    )
    if status != "converged_by_solver_report":
        print(f"CFD convergence: {status}. Finite output may be used for adaptation.", flush=True)
    return {"actual_flow_iterations": actual, "cfd_convergence_status": status,
            "final_history_values": json.dumps(last, sort_keys=True)}


def verify_cfd_outputs(mesh_path, directory, sensor, flow_iterations):
    """Validate the selected field, all restart numbers and mesh association."""
    restart, vtu, history = (directory / name for name in
                             ("restart_flow.csv", "flow.vtu", "history.csv"))
    for path in (restart, vtu, history):
        require_output(path)
    field = resolve_sensor_name(restart, sensor)
    points, triangles, _, _ = read_su2_mesh(mesh_path)
    point_ids, coordinates, values = read_restart_field(restart, field)
    match_scalar_field(points, point_ids, coordinates, values)
    if coordinates is None:
        raise ValueError("Restart coordinates are required to verify the current mesh association.")
    # The elementary extractor checks the selected field. Also refuse a restart
    # with non-finite values in any other physical solution column.
    with restart.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream, skipinitialspace=True)
        next(reader)
        for row in reader:
            while row and not row[-1].strip():
                row.pop()
            if row and any(not math.isfinite(float(value)) for value in row):
                raise ValueError("Non-finite CFD restart data.")
    # PARAVIEW writes appended binary XML. Check its XML prefix without trying
    # to parse the raw appended bytes through an XML-only parser.
    with vtu.open("rb") as stream:
        prefix = stream.read(65536)
    if b'VTKFile type="UnstructuredGrid"' not in prefix:
        raise ValueError("Expected an XML unstructured-grid VTU volume solution.")
    piece = re.search(rb'<Piece\s+NumberOfPoints="(\d+)"\s+NumberOfCells="(\d+)"', prefix)
    if (piece is None or int(piece.group(1)) != len(points)
            or int(piece.group(2)) != len(triangles)):
        raise ValueError("VTU header does not match the current mesh counts.")
    diagnostics = inspect_history(history, directory / "su2.log", flow_iterations)
    return restart, vtu, diagnostics


def determine_plot_limits(mesh_data, chord):
    """Fix all later plot limits from baseline coordinates and airfoil edges."""
    points, _, edges, tags = mesh_data
    minimum, maximum = points.min(axis=0), points.max(axis=0)
    padding = 0.02 * (maximum - minimum)
    airfoil = points[np.unique(edges[tags == "airfoil"])]
    near_min, near_max = airfoil.min(axis=0), airfoil.max(axis=0)
    return (
        (minimum - padding, maximum + padding),
        (near_min - np.array([0.25, 0.5]) * chord,
         near_max + np.array([1.0, 0.5]) * chord),
    )


def plot_mesh(mesh_path, output_path, *, limits, solver, case_name, cycle,
              subiteration, global_step, target_nodes, sensor):
    """Plot actual connectivity and markers with identical limits every time."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    points, triangles, edges, tags = read_su2_mesh(mesh_path)
    figure, axes = plt.subplots(1, 2, figsize=(15, 6), layout="constrained")
    try:
        colours = {"airfoil": "#b94732", "farfield": "#31788e", "outlet": "#7a5598"}
        for axis, title, (minimum, maximum) in zip(
                axes, ("Complete computational domain", "Airfoil and near wake"), limits):
            # triplot uses the supplied SU2 triangles, never a new triangulation.
            axis.triplot(points[:, 0], points[:, 1], triangles,
                         color="#536677", linewidth=0.18, rasterized=True)
            for tag in np.unique(tags):
                axis.add_collection(LineCollection(
                    points[edges[tags == tag]], colors=colours.get(tag, "#555555"),
                    linewidths=0.8, label=tag,
                ))
            axis.set(xlim=(minimum[0], maximum[0]), ylim=(minimum[1], maximum[1]),
                     title=title, xlabel="x [mesh coordinate units]", ylabel="y [mesh coordinate units]")
            axis.set_aspect("equal", adjustable="box")
            axis.legend(loc="upper right", fontsize=8)
        budget = "not applicable" if target_nodes is None else str(target_nodes)
        figure.suptitle(
            f"{solver} | {case_name} | Cycle {cycle}, subiteration {subiteration}, step {global_step}\n"
            f"Vertices: {len(points)} | Triangles: {len(triangles)} | "
            f"Indicative budget: {budget} | Sensor: {sensor}"
        )
        figure.savefig(output_path, dpi=200)
    finally:
        plt.close(figure)


def write_pvd_collection(output_path, datasets):
    """Atomically replace one portable collection of completed solve outputs."""
    root = ET.Element("VTKFile", type="Collection", version="0.1", byte_order="LittleEndian")
    collection = ET.SubElement(root, "Collection")
    for step, path in datasets:
        relative = Path(os.path.relpath(path, output_path.parent)).as_posix()
        ET.SubElement(collection, "DataSet", timestep=str(step), group="", part="0", file=relative)
    ET.indent(root)
    temporary = output_path.with_suffix(".pvd.tmp")
    with temporary.open("wb") as stream:
        ET.ElementTree(root).write(stream, encoding="utf-8", xml_declaration=True)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, output_path)


def append_history(history_path, record):
    """Append only a completed step; missing/unavailable diagnostics stay blank."""
    new_file = not history_path.exists()
    with history_path.open("a", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=HISTORY_FIELDS)
        if new_file:
            writer.writeheader()
        writer.writerow(record)
        stream.flush()
        os.fsync(stream.fileno())


def sha256_file(path):
    """Hash preserved inputs and program sources without loading them at once."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_manifest(run_directory, inputs, dependencies, container_info):
    """Retain provenance; observed image identity is not a solver build proof."""
    sources = [Path(__file__)] + [PROJECT_ROOT / "02_mesh_adaptation" / name for name in
                                ("extract_solution.py", "mesh_adaptation.py", "run_simulation.py")]
    manifest = {
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "python_version": sys.version, "numpy_version": np.__version__,
        "expected_su2_version": "8.5.0", "container_name": CONTAINER_NAME,
        "container_image_id": container_info.get("Image"),
        "container_was_running_at_preflight": container_info["State"]["Running"],
        "input_sha256": {name: sha256_file(path) for name, path in inputs.items()},
        "dependency_sha256": {key: sha256_file(path) for key, path in dependencies.items()},
        "source_sha256": {path.relative_to(PROJECT_ROOT).as_posix(): sha256_file(path) for path in sources},
        "launch_policy": "Synchronous serial SU2_CFD through run_simulation.py; no solution transfer",
        "evidence_status": "execution_record_only; numerical_and_physical_validation_not_established",
    }
    for key, arguments in (("git_commit", ["rev-parse", "HEAD"]),
                           ("git_status", ["status", "--short"])):
        result = subprocess.run(["git", "-C", str(PROJECT_ROOT), *arguments],
                                capture_output=True, text=True, check=True)
        manifest[key] = result.stdout.strip()
    (run_directory / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


def print_workflow_summary(mesh_path, config_path, adaptation_path, parameters,
                            solver, case_directory, run_directory):
    """Show the complete schedule before automatically starting the baseline."""
    print("\n" + "=" * 60 + "\nWORKFLOW CONFIGURATION\n" + "=" * 60)
    for label, value in (("Initial mesh", mesh_path), ("CFD configuration", config_path),
                         ("Adaptation settings", adaptation_path), ("CFD solver", solver),
                         ("Simulation case", case_directory.name), ("Run identifier", run_directory.name),
                         ("Sensor", parameters["SENSOR"]), ("Reference chord", parameters["CHORD"]),
                         ("Number of cycles", len(parameters["SUBITER"])),
                         ("Total adaptations", sum(parameters["SUBITER"])),
                         ("Total CFD simulations", 1 + sum(parameters["SUBITER"]))):
        print(f"{label + ':':24s}{value}")
    for cycle in range(len(parameters["SUBITER"])):
        print(f"\nCycle {cycle + 1}:")
        for key in CYCLE_KEYS:
            print(f"    {key}: {parameters[key][cycle]}")
    for key in ("HMAX", "EPSILON", "COLLAPSE_THRESHOLD", "SPLIT_THRESHOLD", "REMESH_MAX_ITER"):
        print(f"{key}: {parameters[key]}")
    print(f"\nResults directory:\n    {run_directory}\n" + "=" * 60, flush=True)


# ---------------------------------------------------------------------------
# Sequential baseline and two explicit adaptation loops.
# ---------------------------------------------------------------------------
def adaptive_workflow(mesh_path, config_path, adaptation_config_path):
    """Run the complete workflow without prompting; return its run directory.

    All original inputs remain untouched. Each invocation reserves fresh run
    storage; there is no automatic resume and no deletion of partial results.
    """
    started = time.perf_counter()
    mesh_path = resolve_input_file(mesh_path, {".su2"})
    mesh_data = validate_mesh(mesh_path)
    config_path = resolve_input_file(config_path, {".cfg"})
    solver = detect_cfd_solver(config_path)
    adaptation_path = resolve_input_file(adaptation_config_path, {".cfg", ".txt"})
    parameters = read_adaptation_config(adaptation_path)
    dependencies = resolve_dependencies(config_path)

    # Read-only Docker preflight is delegated to the runner's existing helpers.
    # At future runtime this diagnoses Docker Desktop/mount access before CFD,
    # extraction or remeshing; lifecycle changes remain inside run_simulation.
    try:
        container_info = inspect_container(CONTAINER_NAME)
        category = PROJECT_ROOT / ("04_euler_simulations" if solver == "EULER"
                                   else "05_rans_simulations")
        host_to_container_path(category, container_info)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        raise RuntimeError(f"Docker preflight failed; require an accessible writable bind mount: {error}") from error
    check_metric_feasibility(mesh_data, parameters)
    limits = determine_plot_limits(mesh_data, parameters["CHORD"])
    case = resolve_case_directory(PROJECT_ROOT, solver, mesh_path, config_path,
                                  adaptation_path, dependencies)
    # Existing case symlinks must also map to a writable bind mount.
    host_to_container_path(case, container_info)
    run = create_run_directory(case)
    location = "input snapshot preparation"
    try:
        host_to_container_path(run, container_info)
        # Preserve originals separately from the local execution configuration.
        shutil.copy2(config_path, run / "reference_config.cfg")
        shutil.copy2(adaptation_path, run / "adaptation.cfg")
        dependency_paths = snapshot_dependencies(dependencies, run)
        retained_dependencies = {key: run / value for key, value in dependency_paths.items()}
        write_manifest(run, {"mesh": mesh_path, "reference_config": config_path,
                             "adaptation_config": adaptation_path}, retained_dependencies, container_info)
        print_workflow_summary(mesh_path, config_path, adaptation_path, parameters, solver, case, run)
        reference = run / "reference_config.cfg"
        sensor = parameters["SENSOR"]
        collection = []
        total_adaptation_time = 0.0
        total_cfd_time = 0.0
        global_step = 0
        context = {"solver": solver, "case_name": case.name, "run_id": run.name, "sensor": sensor}

        # Baseline is solve zero and uses the first cycle's iteration limit.
        location = "baseline (cycle 0, subiteration 0)"
        step_started = time.perf_counter()
        baseline = run / "baseline"
        baseline.mkdir()
        current_mesh = baseline / "mesh.su2"
        print("\nBASELINE SIMULATION\n" + "-" * 60, flush=True)
        print(f"Initial vertices: {len(mesh_data[0])}\nInitial triangles: {len(mesh_data[1])}\n"
              f"CFD iteration limit: {parameters['FLOW_ITER'][0]}", flush=True)
        print("[1/3] Preparing baseline files", flush=True)
        shutil.copy2(mesh_path, current_mesh)
        local_config = prepare_cfd_config(reference, current_mesh, parameters["FLOW_ITER"][0],
                                          baseline, solver, sensor, retained_dependencies)
        print("[2/3] Generating baseline mesh plot", flush=True)
        timer = time.perf_counter()
        plot_mesh(current_mesh, baseline / "mesh.png", limits=limits, solver=solver,
                  case_name=case.name, cycle=0, subiteration=0, global_step=0,
                  target_nodes=None, sensor=sensor)
        plot_time = time.perf_counter() - timer
        print("[3/3] Running CFD simulation", flush=True)
        timer = time.perf_counter()
        run_simulation(current_mesh, local_config, log_path=baseline / "su2.log")
        cfd_time = time.perf_counter() - timer
        current_solution, vtu, cfd_diagnostics = verify_cfd_outputs(
            current_mesh, baseline, sensor, parameters["FLOW_ITER"][0])
        total_cfd_time += cfd_time
        collection.append((0, vtu))
        append_history(run / "adaptation_history.csv", {
            **context, **cfd_diagnostics, "global_step": 0, "cycle": 0, "subiteration": 0,
            "actual_nodes": len(mesh_data[0]), "actual_triangles": len(mesh_data[1]),
            "flow_iteration_limit": parameters["FLOW_ITER"][0], "plot_runtime_seconds": plot_time,
            "cfd_runtime_seconds": cfd_time, "total_runtime_seconds": time.perf_counter() - step_started,
            "mesh_path": current_mesh.relative_to(run).as_posix(),
            "restart_path": current_solution.relative_to(run).as_posix(),
            "vtu_path": vtu.relative_to(run).as_posix(),
        })
        write_pvd_collection(run / "adaptation.pvd", collection)
        del mesh_data

        # Invariant at each loop entry: current_solution belongs to current_mesh.
        # Only replace this pair after the new solve passes output validation.
        for c in range(len(parameters["SUBITER"])):
            cycle = c + 1
            cycle_directory = run / f"cycle_{cycle:02d}"
            cycle_directory.mkdir()
            cycle_collection = []
            budget = parameters["BUDGET"][c]
            h_min = parameters["HMINS"][c] * parameters["CHORD"]
            h_max = parameters["HMAX"] * parameters["CHORD"]
            for subiteration in range(1, parameters["SUBITER"][c] + 1):
                global_step += 1
                location = f"cycle {cycle}, subiteration {subiteration}, global step {global_step}"
                directory = cycle_directory / f"subiter_{subiteration:02d}"
                directory.mkdir()
                step_started = time.perf_counter()
                print("\n" + "=" * 60 + f"\nCYCLE {cycle}/{len(parameters['SUBITER'])} - "
                      f"SUBITERATION {subiteration}/{parameters['SUBITER'][c]}\n" + "=" * 60)
                print(f"Global adaptation step: {global_step}/{sum(parameters['SUBITER'])}\n"
                      f"Indicative vertex budget: {budget}\nTarget complexity: {target_complexity(budget):.12g}\n"
                      f"Lp exponent: {parameters['P_NORMS'][c]}\nMinimum metric size: {h_min:g}\n"
                      f"Maximum metric size: {h_max:g}\nMaximum aspect ratio: {parameters['ARMAXS'][c]}", flush=True)
                print("[1/4] Extracting CFD solution", flush=True)
                timer = time.perf_counter()
                field = resolve_sensor_name(current_solution, sensor)
                npz_path = directory / "adaptation_input.npz"
                extract_solution(current_mesh, current_solution, field, npz_path)
                extraction_time = time.perf_counter() - timer

                print("[2/4] Adapting mesh", flush=True)
                timer = time.perf_counter()
                new_mesh = directory / "mesh.su2"
                diagnostics = {}
                status = mesh_adaptation(
                    input_path=npz_path, target_nodes=budget, epsilon=parameters["EPSILON"],
                    p=parameters["P_NORMS"][c], h_min=h_min, h_max=h_max,
                    max_aspect_ratio=parameters["ARMAXS"][c],
                    collapse_threshold=parameters["COLLAPSE_THRESHOLD"],
                    split_threshold=parameters["SPLIT_THRESHOLD"],
                    max_iterations=parameters["REMESH_MAX_ITER"], output_path=new_mesh,
                    diagnostics=diagnostics,
                )
                adaptation_time = time.perf_counter() - timer
                if status not in ("converged", "stationary", "iteration_limit"):
                    raise ValueError(f"Unknown remeshing status: {status}")
                adapted_data = validate_mesh(new_mesh)
                nodes, triangles = len(adapted_data[0]), len(adapted_data[1])
                del adapted_data
                if status != "converged":
                    print(f"Remeshing status: {status}; metric remeshing is not fully converged.", flush=True)
                print("[3/4] Generating mesh plot", flush=True)
                timer = time.perf_counter()
                plot_mesh(new_mesh, directory / "mesh.png", limits=limits, solver=solver,
                          case_name=case.name, cycle=cycle, subiteration=subiteration,
                          global_step=global_step, target_nodes=budget, sensor=sensor)
                plot_time = time.perf_counter() - timer
                local_config = prepare_cfd_config(reference, new_mesh, parameters["FLOW_ITER"][c],
                                                  directory, solver, sensor, retained_dependencies)
                print("[4/4] Running CFD simulation", flush=True)
                timer = time.perf_counter()
                run_simulation(new_mesh, local_config, log_path=directory / "su2.log")
                cfd_time = time.perf_counter() - timer
                new_solution, vtu, cfd_diagnostics = verify_cfd_outputs(
                    new_mesh, directory, sensor, parameters["FLOW_ITER"][c])
                step_time = time.perf_counter() - step_started

                # The NPZ describes the previous pair, while mesh.su2 and the
                # outputs here describe the newly solved pair. Record only now.
                record = {
                    **context, **diagnostics, **cfd_diagnostics,
                    "global_step": global_step, "cycle": cycle, "subiteration": subiteration,
                    "target_nodes": budget, "actual_nodes": nodes, "actual_triangles": triangles,
                    "target_complexity": target_complexity(budget), "h_min": h_min, "h_max": h_max,
                    "p": parameters["P_NORMS"][c], "epsilon": parameters["EPSILON"],
                    "max_aspect_ratio": parameters["ARMAXS"][c],
                    "collapse_threshold": parameters["COLLAPSE_THRESHOLD"],
                    "split_threshold": parameters["SPLIT_THRESHOLD"],
                    "remesh_max_iterations": parameters["REMESH_MAX_ITER"],
                    "flow_iteration_limit": parameters["FLOW_ITER"][c], "remeshing_status": status,
                    "extraction_runtime_seconds": extraction_time,
                    "adaptation_runtime_seconds": adaptation_time, "plot_runtime_seconds": plot_time,
                    "cfd_runtime_seconds": cfd_time, "total_runtime_seconds": step_time,
                    "mesh_path": new_mesh.relative_to(run).as_posix(),
                    "restart_path": new_solution.relative_to(run).as_posix(),
                    "vtu_path": vtu.relative_to(run).as_posix(),
                }
                append_history(run / "adaptation_history.csv", record)
                collection.append((global_step, vtu))
                cycle_collection.append((global_step, vtu))
                write_pvd_collection(run / "adaptation.pvd", collection)
                write_pvd_collection(cycle_directory / f"cycle_{cycle:02d}.pvd", cycle_collection)
                current_mesh, current_solution = new_mesh, new_solution
                total_adaptation_time += adaptation_time
                total_cfd_time += cfd_time
                print(f"Actual vertices: {nodes}\nActual triangles: {triangles}\nRemeshing status: {status}\n"
                      f"Adaptation runtime: {adaptation_time:.3f} s\nCFD runtime: {cfd_time:.3f} s\n"
                      f"Total step runtime: {step_time:.3f} s\nStep completed.", flush=True)

        print("\n" + "=" * 60 + "\nWORKFLOW COMPLETED\n" + "=" * 60)
        print(f"CFD solver: {solver}\nSimulation case: {case.name}\nRun identifier: {run.name}\n"
              f"Completed cycles: {len(parameters['SUBITER'])}\nCompleted adaptations: {global_step}\n"
              f"Completed CFD solves: {global_step + 1}\nFinal vertices: {nodes}\nFinal triangles: {triangles}\n"
              f"Total adaptation time: {total_adaptation_time:.3f} s\nTotal CFD time: {total_cfd_time:.3f} s\n"
              f"Total workflow time: {time.perf_counter() - started:.3f} s\n"
              f"Results directory: {run}\nFinal mesh: {current_mesh}\nFinal restart: {current_solution}\n"
              f"Final CFD solution: {vtu}\nGlobal ParaView animation: {run / 'adaptation.pvd'}\n"
              f"Workflow history: {run / 'adaptation_history.csv'}", flush=True)
        return run
    except BaseException as error:
        # KeyboardInterrupt is deliberately included; retain every partial file.
        # Cleanup of a solver/container is exclusively the runner's concern.
        message = f"Workflow stopped at {location}: {type(error).__name__}: {error}\nResults retained: {run}\n"
        print(message, file=sys.stderr, flush=True)
        try:
            (run / "failure.txt").write_text(message, encoding="utf-8")
        except OSError as report_error:
            print(f"Could not retain failure report: {report_error}", file=sys.stderr)
        raise


def main():
    """Expose the three sequential prompts or accept exactly three arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", metavar="PATH",
                        help="Initial mesh, reference CFD config and adaptation config, in that order.")
    arguments = parser.parse_args()
    if len(arguments.paths) not in (0, 3):
        parser.error("Provide all three positional paths or none of them.")
    try:
        if arguments.paths:
            mesh_path, config_path, adaptation_path = arguments.paths
        else:
            print("=" * 60 + "\nAUTOMATED ANISOTROPIC CFD MESH ADAPTATION\n" + "=" * 60)
            mesh_path, mesh_data = request_input_file(
                "\nEnter the initial SU2 mesh path:", {".su2"}, validate_mesh)
            print(f"Mesh successfully loaded.\n\nVertices: {len(mesh_data[0])}\nTriangles: {len(mesh_data[1])}")
            del mesh_data
            config_path, solver = request_input_file(
                "\nEnter the reference SU2 CFD configuration path:", {".cfg"}, detect_cfd_solver)
            print(f"CFD configuration successfully loaded.\n\nDetected solver: {solver}")
            adaptation_path, _ = request_input_file(
                "\nEnter the adaptation configuration path:", {".cfg", ".txt"}, read_adaptation_config)
            print("Adaptation configuration successfully loaded.", flush=True)
        adaptive_workflow(mesh_path, config_path, adaptation_path)
    except (KeyboardInterrupt, EOFError):
        print("\nWorkflow interrupted. Existing results have been retained.", file=sys.stderr)
        return 130
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr, flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
