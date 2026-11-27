# Triangular anisotropic mesh adaptation

Mission 02 implements one adaptation of a planar linear-triangle mesh from nodal density. It follows the [operational equations and verification map](mathematical-specification.md) extracted from the project's [LaTeX source](Theory_Latex/main.tex). See the [verification report](verification.md) for the manufactured cases and the Mission 01 mesh run.

## Reproduce the verification

Run from the repository root using the existing project environment. The numerical core needs Python 3.10 or later and NumPy 1.26 or later. Matplotlib 3.8 or later is required only for plots. The cloud verification used Python 3.12, NumPy 2.5.5 and Matplotlib 3.11.2. The host environment inventory records the user's Python 3.14 setup.

```sh
project_python() {
  .PROJECT_CFDvenv/bin/python -B -c 'import pathlib,sys; p=pathlib.Path(".PROJECT_CFDvenv").absolute(); assert pathlib.Path(sys.prefix)==p and pathlib.Path(sys.executable).parent==p/"bin"' &&
  .PROJECT_CFDvenv/bin/python -B "$@"
}
export PYTHONPATH=src
export OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR=local/mpl-cache
mkdir -p "$MPLCONFIGDIR"
project_python -m unittest discover -s tests -v
project_python 02_mesh_adaptation/tools/verify_adaptation.py \
  outputs/mission-02/verification --diagnostics full --figures
```

The driver requires a new output directory. It writes `results.json`, a summary and an adapted mesh per case, and optional VTU, NPZ and figure outputs. The mesh files and graphics are reproducible run artifacts under the ignored `outputs/` directory. Use `--case square-19 --passes 4 --profile` with a fresh directory for profiling a shorter run. The core uses NumPy. The standard library handles input, output, indexing orchestration and profiling. Neither Gmsh nor a CFD solver supplies any adaptation operation.

## Adapt a mesh and solution

Use the Mission 01 SU2 snapshot or another planar unpartitioned SU2 ASCII mesh. The reader accepts type-5 triangles and type-3 line markers, validates mesh connectivity and reorients clockwise triangles. It preserves explicit SU2 point identifiers for field alignment. The writer uses contiguous point and triangle IDs, with marker tags preserved. Mission 01's Gmsh MSH file records the original geometry entities. The companion SU2 file is the input for this adapter.

The nodal CSV requires a `Density` column plus `PointID`, planar `x` and `y`, or both forms of identity. The reader accepts shuffled rows and checks for missing, duplicate, ambiguous and inconsistent matches. Supply both IDs and coordinates wherever available. By default the geometric match allows an absolute componentwise difference of `1e-9` times the largest mesh extent. A `solution` JSON object can override `density_column`, `id_column`, `coordinate_columns` and `coordinate_atol`. The supported path is an ASCII CSV with one row per source vertex. Native binary restart, partitioned meshes and surface-only solution files require a separate input step.

A minimal `run.json` for the unit-square verification scale is:

```json
{
  "metric": {
    "target_nodes": 280,
    "epsilon": 0.05,
    "h_min": 0.023,
    "h_max": 0.35,
    "max_aspect": 7,
    "p": 1
  },
  "remeshing": {"max_passes": 8},
  "diagnostics": "vtk"
}
```

```sh
project_python -m cfd_adapt.cli source.su2 solution.csv run.json outputs/mission-02/my-run
```

Use a new output directory to retain the provenance of earlier runs. The command writes `adapted.su2`, `summary.json` and a copy of the input configuration. All quality, length and geometry statistics are reported even when termination reaches a budget.

## Parameters

| Input | Default and admissible values | Units and effect |
| --- | --- | --- |
| `target_nodes`, $N_\star$ | Positive integer, choose exactly one budget | Vertices. The formula $C_\star\simeq\sqrt3 N_\star/2$ is approximate. |
| `target_complexity`, $C_\star$ | Positive finite scalar, choose exactly one budget | Dimensionless metric area. Must fall in the feasible interval. |
| `p` | Positive integer, default 1 | Dimensionless $L^p$ objective exponent. Larger values accentuate high-error regions. |
| `epsilon`, $\varepsilon$ | Required, strictly positive finite | Density per length squared. Floor applied to absolute Hessian eigenvalues. |
| `h_min`, `h_max` | Required, $0<h_{\min}\le h_{\max}<\infty$ | Same units as coordinates. Limit the principal sizes of the metric. |
| `max_aspect`, $A_{\max}$ | Required, finite and $\ge 1$ | Dimensionless maximum principal-size ratio. One imposes an isotropic metric. |
| `complexity_rtol` | $10^{-10}$, strictly between zero and one | Dimensionless. Absolute log-complexity error used in bisection. |
| `bisection_steps` | 100, positive integer | Upper bound on log-space complexity evaluations after bracketing. |
| `collapse_below`, $\ell_c$ | 0.62, strictly between zero and one | Dimensionless. Shorter metric edges can be collapsed. |
| `split_above`, $\ell_s$ | 1.45, greater than one and finite | Dimensionless. Longer metric edges can be split. |
| `min_quality` | 0.08, positive and at most `convergence_quality` | Local metric quality guard for split, collapse and motion. If the current patch already falls below it, the guard uses the current patch minimum. |
| `convergence_quality` | 0.20, at most one | Final minimum metric mean ratio for reporting `criteria_met`. |
| `quality_gain` | $10^{-5}$, strictly between zero and one | Dimensionless local improvement threshold and relative allowance in patch edge checks. |
| `relocation_steps` | 4, positive integer | At most this many backtracking trial positions per vertex per pass. |
| `max_passes`, `max_operations` | 8 and 100,000, positive integers | Sweep and accepted-operation budgets. A full sweep scans local candidates once per phase. |
| `location_tolerance` | $10^{-10}$, strictly between zero and $10^{-4}$ | Dimensionless barycentric and segment-parameter allowance for source-boundary round-off. |
| `diagnostics` | `summary`, `vtk` or `full` | Adds source/adapted triangle and edge VTUs at `vtk` and compressed raw arrays at `full`. |

Remeshing defaults are implementation choices for the verification cases. They require assessment for CFD workloads. Floating-point degeneracy screening uses a signed double-area threshold of 64 machine epsilons times the squared longest local edge. All lengths use a consistent input coordinate scale. Raw Hessians have units of density per length squared. Metric tensors have units of inverse length squared. The prescribed size limits govern the field eigenvalues. The resulting mesh can still contain short and long physical edges.

## Local operations and geometry

- A sweep tries the shortest edges first for collapse, the longest edges first for split, then swaps and interior vertex relocation. Edge owners, vertex fans, neighbors and length caches update with each edit. Candidate lengths are rechecked before acceptance.
- A split bisects the integrated metric length and replaces all incident triangles on the edge when the children pass orientation and metric quality checks. A collapse checks the vertex-link condition, positive patch areas, local metric quality, and newly introduced long edges. It removes only interior, unprotected vertices. It ranks the two endpoint choices by local metric quality and preserves the coordinates of the surviving endpoint.
- A swap requires a convex patch, a new diagonal and a strict improvement in minimum metric mean ratio. The new diagonal must meet the length band or reduce the old diagonal's deviation on its side of the band. Relocation proposes a neighbor-weighted target using the sampled metric tensors. It backtracks, requiring positive areas, a local shape guard, and improvement in a score combining minimum and mean metric quality with mean squared log metric length. It checks candidate edge lengths against the same immutable source.
- For triangle quality, $q=4\sqrt3 A/(a^2+b^2+c^2)$ in physical space. Metric mean ratio applies the same expression with a constant tensor equal to the mean of the three source-field tensors sampled at the triangle vertices. Edge threshold checks use the exact piecewise-P1 integral. The two diagnostics measure different properties. An equilateral triangle has mean ratio one.
- A new vertex cannot be removed in the next collapse phase. Edge length bands, quality improvement thresholds and operation budgets limit reversals. `criteria_met` requires a stationary sweep, all metric edge lengths in range, and the specified final minimum metric quality. `stationary` indicates remaining violations after all eligible operations were rejected. `pass_limit` and `operation_limit` record exhaustion of budgets. The summary reports short and long edges, constrained violations, rejected operations by reason and mesh changes per pass.

The default geometry policy keeps all boundary vertices and segments fixed. Users can allow splits on verified straight CAD segments through `constraints.splittable_point_id_edges` in the CLI, or permit splitting the entire discrete boundary with JSON `null` for a polygonal geometry. A boundary split lies on the input segment, inherits its marker and keeps both endpoints. For curved CAD entities keep their source segments fixed. Projecting split points onto a CAD curve can take them outside the polygon where the source field is defined, so this policy preserves the existing CAD samples and the original polygonal approximation. A future CAD reprojection feature would require an explicit metric extension or geometry-aware mapping.

Use `constraints.fixed_point_ids` to pin existing features and `constraints.protected_point_id_edges` to retain feature edges as well as their endpoints. Python `GeometryConstraints` accepts source array indices, while the CLI accepts external point IDs. For the Mission 01 evidence the driver identifies straight CAD segments, keeps the curved airfoil and outer arcs fixed, and protects the embedded wake edges. Those edges remain interior. The reader and the adaptation core apply no airfoil-specific logic. The distinction matters when running from SU2, because its physical markers do not encode the embedded wake.

## Python interface and Mission 03

```python
from cfd_adapt import GeometryConstraints, MetricParameters, RemeshParameters, adapt
from cfd_adapt.io import read_density_csv, read_su2, write_su2

source = read_su2("source.su2")
density = read_density_csv("solution.csv", source)
controls = MetricParameters(epsilon=0.05, h_min=0.023, h_max=0.35,
                            max_aspect=7, target_nodes=280)
result = adapt(source, density, controls, RemeshParameters(), GeometryConstraints())
write_su2("adapted.su2", result.mesh)
```

The API returns the adapted mesh, input density, P1 gradients, raw and symmetric Hessians, all metric tensors and spectra, per-edge and per-triangle diagnostics, timings, and source vertex indices. Retained indices may have moved. New split vertices map to `-1`. Output numbering is compact. An index mapping alone is insufficient for transferring a solution to relocated or new vertices. Mission 03 will own any transfer or restart policy, solver invocation and comparison of density at common physical sample points.
