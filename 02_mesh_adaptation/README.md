# Elementary extraction, adaptation and CFD programs

These three authoritative programs are called directly by [adaptive_workflow.py](../03_adaptive_workflow/adaptive_workflow.py). They remain in this directory; no package relocation or duplicate implementation is required. Their presence does not establish runtime or numerical validation. Mission 03 made only the compatibility changes described below and performed static inspection.

## Solution extraction

```python
extract_solution(mesh_path, restart_path, field_name, output_path)
read_su2_mesh(mesh_path)
read_restart_field(restart_path, field_name)
match_scalar_field(points, point_ids, coordinates, values)
```

[extract_solution.py](extract_solution.py) requires a 2D ASCII SU2 mesh, an ASCII CSV restart, the exact scalar column name and an `.npz` destination. It writes the seven arrays below and returns `None`. Its reader returns `(points, triangles, boundary_edges, boundary_tags)`; only volume type 5 triangles and boundary type 3 edges are supported. Connectivity indexes the coordinate array directly with contiguous zero-based IDs. Original numbering/order and coordinates are preserved. Quadrilaterals, partitioned meshes and 3D meshes are unsupported.

The restart must contain `PointID` and the selected scalar. `x` and `y` must both be present or both absent. The elementary matcher checks complete unique IDs, finite values and coordinate agreement when available; the workflow additionally requires coordinates to substantiate mesh association. It resolves sensor spelling case-insensitively before passing the exact header name. Missing quantities are never reconstructed. Only NumPy and the standard library are required. The extraction source was not changed in Mission 03.

| NPZ array | Shape | Meaning |
| --- | --- | --- |
| `points` | `(N, 2)` | Coordinates in mesh coordinate units |
| `triangles` | `(T, 3)` | Vertex indices into `points` |
| `scalar` | `(N,)` | Scalar value at the same indexed vertex |
| `boundary_edges` | `(E, 2)` | Physical boundary connectivity |
| `boundary_tags` | `(E,)` | Marker name for each boundary edge |
| `boundary_nodes` | `(P,)` | Unique protected vertex indices |
| `field_name` | scalar string | Exact extracted column name |

No scalar/coordinate normalization or solution interpolation is performed. NPZ is loaded without pickle by the adaptation program.

## Triangular mesh adaptation

```python
mesh_adaptation(input_path, target_nodes, epsilon, p, h_min, h_max,
                max_aspect_ratio, collapse_threshold, split_threshold,
                max_iterations, output_path='adapted_mesh.su2', *,
                diagnostics=None)
```

[mesh_adaptation.py](mesh_adaptation.py) consumes the seven-array NPZ and explicit metric/remeshing parameters, writes a new SU2 mesh and returns `converged`, `stationary` or `iteration_limit`. `p`, `target_nodes` and `max_iterations` are positive integers; lengths use mesh coordinate units, epsilon uses scalar-Hessian units, and edge thresholds are dimensionless metric lengths. Dependencies are NumPy and the standard library.

The existing method uses double P1 recovery, Hessian symmetrisation, absolute eigenvalues with an epsilon floor, the analytical Lp metric, size/aspect-ratio bounding and complexity restoration. Its continuous P1 source metric remains fixed throughout one remesh. Local operations split/collapse edges, swap diagonals and relocate vertices. The workflow recomputes this metric from the latest CFD field at each adaptation step; it adds no alternative reconstruction or remesher.

`target_nodes` is an indicative vertex budget: the actual mathematical target is `C_star = sqrt(3)/2 * target_nodes`. Bounding restores metric complexity; it neither imposes the final vertex count nor proves remeshing/spatial convergence. A stationary or iteration-limited valid mesh can be used, with its status explicitly retained.

The optional `diagnostics` dictionary is a backward-compatible observation interface added in Mission 03. The internal `adapt_mesh(..., max_iterations, diagnostics=None)` forwards already computed final metric/physical quality statistics, metric-edge RMS deviation from one, short/long edge counts and frozen-edge violations; operation counts accumulate across passes. The public wrapper adds bounded complexity and relative complexity error. Neither return contract nor numerical choices change. No value is inferred from unavailable diagnostics.

### Boundary policy

`airfoil` and `farfield` vertices are protected. Their coordinates and frozen edges are preserved by the existing algorithm. Outlet-only vertices are not protected; outlet edits remain constrained to the existing straight boundary geometry. Shared outlet/farfield endpoints remain protected. The `outlet` marker and all physical boundary markers remain necessary in the exported mesh. The algorithm also freezes non-outlet boundary edges, checks boundary topology and area conservation, and rejects invalid final triangles. The workflow neither rebuilds the boundary nor changes this policy.

No hybrid boundary-layer or quadrilateral support is implemented. Large-mesh performance, mathematical code verification and complete adaptation behavior remain to be assessed separately; inspecting safeguards does not establish their runtime correctness.

## Docker CFD runner

```python
run_simulation(mesh_path, config_path, *, log_path=None)
```

[run_simulation.py](run_simulation.py) requires an existing mesh and `.cfg` accessible through a writable bind mount of `su2-docker`. It returns `None`; SU2 writes into the configuration directory, subject to configured filenames. It uses only the standard library and the Docker CLI. It inspects the existing container, selects the most specific matching host bind mount, rejects read-only access, maps paths, creates a temporary local configuration and runs synchronous serial `SU2_CFD`.

If initially stopped, the container is started for the solve and stopped in `finally`. An already-running container is left running. Temporary launch configuration/PID markers are removed; scientific inputs and outputs are retained. Extraction/adaptation/plotting execute on the host, with no newly started container left running between solves. No image, mount, dependency or resource setting is changed.

Mission 03 adds the documented non-interactive Python/MPI activation, optional terminal-plus-file output streaming, writable-mount checking and safe relative mesh naming for sibling workflow snapshots. On an exception, a unique invocation PID marker enables a guarded SIGTERM attempt against only the matching solver executable and working directory. A shutdown timeout is explicitly reported as unconfirmed; signal propagation and cleanup races require future Docker testing. Container start/stop failures may still require manual diagnosis. No MPI CFD parallelism or runtime timeout is introduced.

`log_path`, when supplied, must be a new file. The original two-argument call and CLI remain accepted. The runner makes no CFD-convergence claim from a zero exit code; the workflow inspects outputs, history and explicit solver messages. Non-interactive activation uses the container-specific path recorded in [environment.md](../docs/environment.md), not the high-level `SU2` Python package.

See the [workflow guide](../03_adaptive_workflow/README.md) for execution instructions and the deferred validation list. No elementary program, Docker command or CFD calculation was executed during Mission 03.
