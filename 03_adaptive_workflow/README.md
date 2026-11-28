# Automated anisotropic CFD mesh adaptation

[adaptive_workflow.py](adaptive_workflow.py) orchestrates the three [elementary programs](../02_mesh_adaptation/README.md). Source implementation and static review are complete; execution, Docker/CFD integration, numerical validation and performance assessment remain pending. Mission 03 permits source/documentation work only. The commands below are for a separate testing and debugging session with the user; they were not executed during implementation.

## Future use

From the checkout root, select the existing project environment and check its identity before executing project code, as required by [environment.md](../docs/environment.md):

```bash
source .PROJECT_CFDvenv/bin/activate
python -B -c 'import pathlib, sys; expected = pathlib.Path(".PROJECT_CFDvenv").absolute(); assert pathlib.Path(sys.prefix) == expected; assert pathlib.Path(sys.executable).parent == expected / "bin"; print(sys.executable, sys.prefix)'
python 03_adaptive_workflow/adaptive_workflow.py
```

The program also checks interpreter identity before importing NumPy. Run this from your actual checkout; no machine-specific absolute project path is required.

`main()` requests and validates these files **sequentially**, retrying an invalid file before the next prompt:

```text
Enter the initial SU2 mesh path:
>
Enter the reference SU2 CFD configuration path:
>
Enter the adaptation configuration path:
>
```

The mesh must be a readable `.su2` with 2D triangles and an `airfoil` boundary marker. After mesh validation, counts are printed. The `.cfg` is parsed and its solver printed before the adaptation `.cfg`/`.txt` is requested. Absolute/relative paths, `~`, spaces, quoted paths and terminal-escaped dragged paths are accepted without shell evaluation. Relative user paths refer to the launch directory; all outputs are located independently from the script's project root.

After all inputs are valid, the future workflow automatically prints its schedule and starts the baseline solve, with no extra solver/output selection or confirmation prompt.

For non-interactive use, provide all three paths or none; partial argument lists are rejected:

```bash
python 03_adaptive_workflow/adaptive_workflow.py \
    /path/to/baseline.su2 \
    /another/path/to/config.cfg \
    /another/path/to/adaptation.cfg
```

These are placeholder paths, not provided cases. No example simulation case or standalone adaptation configuration is included.

The importable interface is:

```python
adaptive_workflow(mesh_path, config_path, adaptation_config_path)
```

It returns the created run directory, never prompts, and propagates errors. Because numbered directories are preserved, callers can load the file with standard `importlib.util.spec_from_file_location`, then call the function using the designated host interpreter. Module import checks the host environment and loads the elementary functions, but starts no solver or interactive input. `PROJECT_ROOT` is derived from `__file__`; the elementary script directory is added explicitly to `sys.path`. No package framework or custom class is introduced.

## Complete adaptation configuration

Every key is required. This example specifies syntax and a proposed schedule; it is not a validated numerical recommendation. In particular, `EPSILON` and `REMESH_MAX_ITER` are provisional starting values.

```ini
% ============================================================
% AUTOMATIC ANISOTROPIC MESH ADAPTATION
% ============================================================

% Number of adaptation steps per cycle
SUBITER = (3, 3, 3, 3)

% Indicative vertex budget per cycle
BUDGET = (40000, 60000, 100000, 150000)

% Maximum CFD solver iterations per cycle
FLOW_ITER = (3000, 3000, 3000, 6000)

% Scalar adaptation sensor
SENSOR = (MACH)

% Reference chord length in mesh coordinate units
CHORD = 1.0

% Maximum relative metric size
HMAX = 100.0

% Minimum relative metric size per cycle
HMINS = (1e-5, 5e-6, 5e-6, 5e-6)

% Lebesgue norm exponent per cycle
P_NORMS = (2.0, 4.0, 8.0, 8.0)

% Maximum metric aspect ratio per cycle
ARMAXS = (1000, 2000, 5000, 500)

% Hessian curvature regularisation floor
EPSILON = 1e-8

% Dimensionless metric edge length thresholds
COLLAPSE_THRESHOLD = 0.62
SPLIT_THRESHOLD = 1.45

% Maximum internal remeshing passes per adaptation step
REMESH_MAX_ITER = 20
```

| Parameter | Meaning and validation |
| --- | --- |
| `SUBITER` | Positive integer adaptation steps per cycle; its nonempty list determines cycle count |
| `BUDGET` | Positive integer indicative vertex budget per cycle; target complexity is `sqrt(3)/2 * BUDGET` |
| `FLOW_ITER` | Positive integer maximum steady `ITER` per solve in that cycle |
| `SENSOR` | Exactly one scalar CSV field name; matching ignores case and rejects ambiguity |
| `CHORD` | Positive reference chord in the mesh coordinate units; coordinates are never rescaled |
| `HMAX` | Positive global maximum metric size relative to chord |
| `HMINS` | Positive per-cycle minimum metric sizes relative to chord; each must be at most `HMAX` |
| `P_NORMS` | Positive integer-valued Lp exponents; `2.0` is accepted as integer 2; fractional exponents are rejected |
| `ARMAXS` | Per-cycle maximum principal metric size ratios, each at least one |
| `EPSILON` | Positive scalar-Hessian curvature floor, passed unchanged |
| `COLLAPSE_THRESHOLD` | Dimensionless metric edge threshold strictly between zero and one |
| `SPLIT_THRESHOLD` | Dimensionless metric edge threshold strictly above one |
| `REMESH_MAX_ITER` | Positive integer maximum internal remeshing passes per adaptation step |

All cycle lists must have equal length. Numbers must be finite; chord-scaled lengths must remain positive/finite. Scalar values, comma-separated parenthesised lists, scientific notation, whitespace, blank lines and `%` full-line/inline comments are supported. Unknown, duplicate, missing or malformed assignments are rejected, without `eval()`. Assignments must fit on one line; multiline continuations and literal `%` inside values/filenames are unsupported. Integer-only controls require integer syntax; only `P_NORMS` additionally accepts integer-valued floating-point syntax.

For cycle `c`, `h_min = HMINS[c] * CHORD` and `h_max = HMAX * CHORD`. Edge thresholds and epsilon are not multiplied by chord. The existing `bounded_metric` feasibility logic is called before CFD using total baseline area and an identity probe; its size-bound interval applies independently of the future scalar field. Extreme numeric magnitudes can still exceed floating-point range and require diagnosis.

An indicative budget is never imposed as an exact final node count. Target complexity, actual node count and available bounded-metric diagnostics are recorded separately. Restored complexity does not establish remeshing convergence.

## Schedule and solver configuration

A cycle holds its parameters constant for `SUBITER[c]` consecutive adaptations. `SUBITER` counts extract/adapt/solve operations; `REMESH_MAX_ITER` limits the internal remeshing passes inside **each** adaptation. The example has four cycles, twelve adaptations and thirteen CFD solves.

1. Copy and plot the baseline; solve with `FLOW_ITER[0]`.
2. For each cycle, and then each subiteration: extract the latest solved field into NPZ, recompute the metric and adapt, read/check the new mesh, plot it, and solve it with `FLOW_ITER[c]`.
3. Verify and record the new mesh/solution pair before using it in the next step, including at a cycle transition. The last adapted mesh receives a final CFD solve.

`SOLVER=EULER` routes to `04_euler_simulations`; `SOLVER=RANS` routes to `05_rans_simulations`. Omitted `SOLVER` defaults to explicitly written Euler in local solve configurations. Commented definitions are ignored; duplicate active assignments and other solvers are rejected. Only steady, fixed-mesh, single-zone direct cases using SU2 mesh input are supported. Unsteady, moving/deforming mesh, adjoint and multizone configurations require a later extension, without changing their physical interpretation silently.

Every local configuration derives from the retained reference. Physical conditions, models, schemes, CFL, boundary conditions, reference quantities and convergence criteria remain unchanged. Execution/output changes set the selected solver, current mesh, maximum `ITER`, fresh initialisation (`RESTART_SOL=NO`, `RESTART_ITER=0`), ASCII restart settings and local output names. No old-mesh restart or solution transfer is used, including for the baseline.

The workflow retains other output formats and adds `RESTART_ASCII` and `PARAVIEW`. If an explicit frequency list exists, it must contain one positive integer for each original format; new formats receive the current iteration limit. Required output groups include `COORDINATES`, `SOLUTION`, `PRIMITIVE`; `COMPACT` is removed and `WRT_RESTART_COMPACT=NO` is written. Existing useful volume/history groups are retained; disabling `NONE` entries are removed, and `ITER` plus every-iteration inner history output enable counting. Auxiliary output paths are assigned distinct local stems to avoid collisions. The intended final files are `restart_flow.csv`, `flow.vtu` and `history.csv`, plus the retained `config.cfg` and streamed `su2.log`.

These settings were checked statically against the [SU2 8.5.0 template](https://github.com/su2code/SU2/blob/v8.5.0/config_template.cfg), [output dispatch](https://github.com/su2code/SU2/blob/v8.5.0/SU2_CFD/src/output/COutput.cpp), [ASCII writer](https://github.com/su2code/SU2/blob/v8.5.0/SU2_CFD/src/output/filewriter/CSU2FileWriter.cpp), [VTU writer](https://github.com/su2code/SU2/blob/v8.5.0/SU2_CFD/src/output/filewriter/CParaviewXMLFileWriter.cpp) and [format map](https://github.com/su2code/SU2/blob/v8.5.0/Common/include/option_structure.hpp). Upstream source review does not establish compatibility with the installed binary.

Mach, Pressure and Temperature belong to `PRIMITIVE`, and Density to `SOLUTION` in [compressible output source](https://github.com/su2code/SU2/blob/v8.5.0/SU2_CFD/src/output/CFlowCompOutput.cpp). Other sensors must be supplied by the reference configuration's appropriate `VOLUME_OUTPUT` groups. CSV display names need not equal SU2 option identifiers; unknown mappings are not invented. Missing/ambiguous fields cause an explicit error listing available columns; no fallback quantity is substituted.

## Cases, snapshots and results

Use the case's own mesh, reference CFD configuration and adaptation configuration. Different source directories and files outside the Docker bind mount are accepted. Relative external file options are resolved against the **original** CFD configuration directory before copying; dependencies receive safe local `inputs/input_NNN.ext` paths in case/run/solve snapshots. Original inputs are never moved or overwritten.

Explicit scalar `FILENAME`/`_FILE` options are handled; restart inputs are unused because restart is disabled. File lists, nested `.cfg` dependencies and explicit runtime configuration reloads are rejected rather than silently changing paths/conditions. Implicit secondary files referenced inside data files or unconventional path-bearing options require review before use; arbitrary SU2 configurations are not claimed supported.

A configuration already under the correct solver category reuses its first-level case directory. Inputs inside an existing `results` directory are rejected. Otherwise the original configuration directory supplies a sanitized case name; directly at a category root, its filename stem is used. Any existing name is preserved and a new suffixed case is reserved. Such new cases receive `baseline.su2`, `config.cfg`, `adaptation.cfg`, the exact `reference_config.cfg` and dependencies. Existing cases are not overwritten when inputs come from separate locations; each run holds the actual selected input snapshot.

Run numbers are independent per case, increase beyond all existing `run_NNN` names and are reserved without overwriting. Future output structure is:

```text
case_name/
├── baseline.su2
├── config.cfg
├── adaptation.cfg
└── results/
    ├── run_001/
    │   ├── reference_config.cfg
    │   ├── adaptation.cfg
    │   ├── inputs/                         # when external files are used
    │   ├── run_manifest.json
    │   ├── adaptation_history.csv
    │   ├── adaptation.pvd
    │   ├── baseline/
    │   │   ├── mesh.su2
    │   │   ├── mesh.png
    │   │   ├── config.cfg
    │   │   ├── su2.log
    │   │   ├── restart_flow.csv
    │   │   ├── flow.vtu
    │   │   └── history.csv
    │   └── cycle_01/
    │       ├── cycle_01.pvd
    │       ├── subiter_01/                 # same solve files, plus adaptation_input.npz
    │       └── subiter_02/
    └── run_002/
```

The NPZ in a subiteration describes its **input** mesh/solution; the sibling `mesh.su2`, restart and VTU describe its **newly solved** mesh. The manifest retains UTC creation time, original-input/dependency/source SHA-256 values, Git commit/status, host Python/NumPy versions and inspected container image identity. It records expected SU2 version, not an asserted runtime build qualification. Physical settings remain in configuration snapshots. Git-ignored results still need an explicit retention decision for a scientific study.

## Visualization and history

Each `mesh.png` uses Matplotlib's non-interactive Agg backend and the actual supplied triangular connectivity. Two equal-scale views show the complete domain/boundaries and the airfoil/near wake. Limits are fixed from baseline coordinates, the `airfoil` marker and `CHORD`; titles report solver, case, cycle, subiteration, global step, mesh counts, indicative budget and sensor. Collection/vectorized drawing avoids an object per triangle; figures are closed after saving. Plot memory/time on large meshes remain unmeasured.

`adaptation.pvd` includes the baseline at index zero and all completed adaptations in order. Each `cycle_NN.pvd` includes only that cycle's completed solves, retaining global step indices. Timesteps are adaptation indices, not physical time. Paths are relative to each collection; VTUs are not duplicated or interpolated. Standard-library XML generation uses atomic replacement after every successful solve, so earlier valid collections survive an interruption. The example schedule ultimately gives 13 global datasets and 3 per cycle.

`adaptation_history.csv` appends one completed-step record, including baseline. It records solver/case/run, indices, sensor, budget and actual counts, target complexity, physical metric bounds, Lp/epsilon/aspect ratio, edge thresholds, remeshing pass/CFD iteration limits and status, extraction/adaptation/plot/CFD/total step timings, and relative mesh/restart/VTU paths. Optional already-computed diagnostics include bounded complexity/error, metric and physical min/median quality, metric-edge RMS, short/long edges, frozen boundary violations and cumulative split/collapse/swap/relocation counts. Baseline-inapplicable entries are blank. Additional columns retain actual inner iterations, reported convergence status and final history values. Times use `perf_counter`; CFD timing includes the runner's lifecycle and output streaming, and total workflow timing includes preparation/verification overhead.

## Failures, evidence and limitations

Fresh solve directories prevent stale output acceptance. Every solve requires readable nonempty restart/VTU/history, a unique selected field, finite restart/history values, complete point IDs and matching restart coordinates. VTU XML prefix and mesh counts are screened without decoding binary payload; complete VTU/ParaView integrity remains unverified.

The [single-zone driver source](https://github.com/su2code/SU2/blob/v8.5.0/SU2_CFD/src/drivers/CSinglezoneDriver.cpp) distinguishes convergence and iteration-limit messages. These are retained as solver reports, separate from physical validation. History's zero-based `Inner_Iter` is converted to a performed iteration count. Unknown outcomes are labelled `unconfirmed`; finite iteration-limited output may proceed with an explicit warning. Detected divergence, non-finite values, reported solver interruption, inconsistent mesh/solution or missing output stops the workflow. A valid `stationary`/`iteration_limit` remesh may proceed with a recorded warning; a budget discrepancy alone is not failure.

The runner owns Docker preflight helpers, activation, execution and lifecycle. The workflow never starts/stops containers itself. A previously running container remains running to protect unrelated work. Failed/interrupted steps retain partial directories and a `failure.txt` when possible; no incomplete solve is registered or recorded as completed. Completed history/PVD entries remain available. There is no automatic resume, parallel solver dispatch, database or GUI. Filesystem failures between CSV/PVD updates may leave collections at different last-completed indices; prior XML remains valid, but cross-file publication is not transactional.

Pending checks include installed SU2 filenames/output groups and final-output behavior, selected scalar availability, non-interactive Docker/mount/signal handling, dependency-path coverage, full VTU readability, large-mesh remeshing/plot performance and end-to-end numerical behavior. No runtime success or convergence claim follows from implementation or static review. RANS routing supports only the same 2D triangular remesher; hybrid/quadrilateral boundary-layer adaptation is a later phase. Runtime testing and debugging must be a separate session with the user.
