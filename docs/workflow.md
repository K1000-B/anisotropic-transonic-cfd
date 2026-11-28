# Working and publication workflow

## Repository organization

The seven numbered stage directories organize scientific responsibilities. The baseline stage contains native geometry, a specialized verifier and a reviewed reference snapshot. Stage 02 contains the three elementary programs; its earlier LaTeX theory is retained in Git history, not the current tree. Stage 03 contains the procedural orchestrator and guide. No simulation case or adaptive output was created in Mission 03. Add cases and tests only in an authorized subsequent phase.

`AGENTS.md` owns permanent collaboration rules; `scientific-framework.md` owns the scientific scope and methodology; `environment.md` owns tooling observations and environment requirements; `project-status.md` owns current progress, decisions and next actions. The README is the public entry point, not a duplicate project log.

### Stage interactions and shared software

`00_resources` supports all stages. `01_baseline_mesh` owns reusable geometry and baseline mesh definitions. `02_mesh_adaptation` owns the adaptation method and its verification evidence. `03_adaptive_workflow` coordinates those tools and SU2. `04_euler_simulations` and `05_rans_simulations` supply case-specific inputs to the same workflow. `06_experimental_validation` compares approved results with traceable experimental evidence. RANS mesh construction may reuse baseline geometry and the Euler-adapted exterior, subject to adequacy checks.

Numbered directories are scientific stages, not Python package names. Mission 03 explicitly preserves `02_mesh_adaptation/extract_solution.py`, `mesh_adaptation.py` and `run_simulation.py` as authoritative sources. `03_adaptive_workflow/adaptive_workflow.py` derives the project root from `__file__`, adds the elementary script directory to `sys.path` and imports their functions directly. There is no `src/` package, custom class or new testing architecture. Do not relocate or duplicate these programs merely for imports. The orchestrator exposes `adaptive_workflow(mesh_path, config_path, adaptation_config_path)` without prompts; interactive input belongs only to `main()`.

Keep dependencies directed: numerical mesh/adaptation code must not import Euler/RANS cases or orchestration; orchestration calls reusable components; cases provide data/configuration; analysis consumes outputs. Share geometry by referencing its authoritative `.geo` definition. ParaView scripts use its application interpreter; host numerical code uses the designated project venv. The solve/extract/adapt pipeline is implemented but unexecuted; no runtime integration claim follows.

### Implemented execution contract (pending execution)

```text
Initial SU2 mesh + Reference CFD configuration + Adaptation configuration
                                |
                                v
                       adaptive_workflow.py
                                |
                                v
                         run_simulation.py
                                |
                                v
                  SU2 CFD resolution inside Docker
                                |
                                v
                 ASCII restart solution + VTU + history
                                |
                                v
                        extract_solution.py
                                |
                                v
                        Seven-array NPZ
                                |
                                v
                        mesh_adaptation.py
                                |
                                v
                        Adapted SU2 mesh
                                |
                                v
                        Mesh visualization
                                |
                                v
                  New CFD resolution inside Docker
                                |
                                v
               Repeat until all adaptation cycles are complete
```

Native macOS Python performs path/config parsing, snapshots, extraction, adaptation, plots, CSV/PVD and provenance. Only CFD resolution executes in the existing container. The runner owns mount discovery, non-interactive activation, synchronous serial execution, logs and container lifecycle; orchestration reuses its read-only preflight helpers without duplicating start/stop logic. Already-running containers remain running; a container started for a solve is stopped before returning to host processing.

The baseline solve uses `FLOW_ITER[0]`. The two explicit loops traverse cycles and their `SUBITER[c]` adaptations; every adaptation uses the latest mesh and its latest verified CFD restart. The new pair replaces the previous pair only after output checks. Cycle transitions preserve that pair and change only the requested cycle controls. The last adapted mesh also receives a CFD solve. No solution transfer is implemented; all solves start from configured reference initialisation with restart disabled. `ITER` is a maximum, not proof of iterative convergence.

Each run records solver-reported convergence separately from remeshing status and scientific validation. A finite iteration-limited solution or valid stationary/iteration-limited remesh may proceed with an explicit record; missing/inconsistent/non-finite outputs stop execution. Collections/history are updated after accepted solves, and partial results are retained on failure/interruption. See the [user guide](../03_adaptive_workflow/README.md) for the complete configuration, output tree, limitations and future commands.

Mission 03 permits static source review only: no Python execution, imports, compilation, tests/fixtures, Docker/SU2 commands, numerical results or dependency changes. Future testing and debugging require a separate session with the user. Preserve the existing LaTeX theory and local abandonment archive.

## Git and GitHub

Use `main` as the primary branch. For now, small verified initialization/documentation changes can go directly to it. Use focused branches and pull requests when changes merit scientific review, involve collaborators or would disrupt the main workflow. Use tags for meaningful reproducible milestones, not every edit.

Before work, inspect status, branch, recent history and remotes. Before synchronization, fetch remote state and inspect divergence; preserve remote contributions. If concurrent changes exist, reconcile non-destructively; request arbitration for substantive scientific conflicts. Never force-push or rewrite history without approval.

Before committing:

1. Inspect every changed file, including untracked content, and the staged diff.
   Before any staging, verify `.PROJECT_CFDvenv/`, `.vscode/` and `local/` are ignored with `git check-ignore`; ensure none are already tracked.
2. Perform checks appropriate to the authorized phase and record outcomes. Mission 03 uses source review, link/path inspection, whitespace and publication review only; no tests/computational programs may be executed. Do not bypass a mandatory hook that would violate this restriction; report the blocker. Later missions may authorize runtime checks separately.
3. Stage explicit paths. Review for sensitive information, redistribution rights and unintended large files. Ignore rules are not a secret scanner.
4. Commit one coherent objective with an English Conventional Commit message (`docs:`, `chore:`, `feat:`, `fix:`, `test:` as appropriate).
5. Push ordinary publishable commits to the configured upstream. Verify local/remote commit equality and final status. If access fails, retain local commits and report the blockage accurately.

The repository is intended to be public. Ordinary commits and pushes are authorized; a license grant, destructive operations and uncertain-rights publication still require user approval. Authentication remains in the existing credential mechanism, never in versioned files. Use a repository-local public GitHub noreply address for commit privacy; do not change global identity settings.

## Existing work and reference material

Before importing external code, course documents or datasets, establish origin, authorship, permitted use and redistribution rights. When uncertain, keep them outside the Git working tree and ask the user. Preserve originals and history when available; record justified modifications separately. Link to protected scientific publications instead of redistributing them. Record bibliographic metadata and persistent identifiers only after verification.

Version public reference identifiers, verified metadata and original summaries in `00_resources/`. For protected originals outside the checkout, use an ignored `local/resource-paths.tsv` mapping from a real reference identifier to an external file path. Create entries only when documents actually arrive. The map is local configuration, not part of the public bibliography or an archival copy of the documents.

## Data and generated results

Version source code, native geometry programs, solver configurations, justified parameters, scripts, tests and publishable documentation. Small, important reference meshes or datasets may be tracked after explicit size, provenance and rights review.

Root `runs/` and `outputs/` remain ignored locations for baseline/other documented computations and derived report figures/tables. Adaptive CFD runs instead belong in `04_euler_simulations/<case>/results/run_NNN/` or the equivalent RANS directory. These case-specific results are ignored; user-authored input meshes/configurations are not blanket-ignored. No generated result belongs in the reusable program directories `02_mesh_adaptation/` or `03_adaptive_workflow/`. Keep inputs and reviewed reference data in their responsible stage, distinguish originals from transformed data, and retain reviewed publishable summaries/manifests with case documentation.

Root `local/` holds ignored machine notes and tool/path configuration; confidential/uncertain-rights originals remain outside the checkout. `.PROJECT_CFDvenv/` and existing `.vscode/` settings are ignored and preserved. No `runs/` or `outputs/` directory is created before real outputs exist. Do not blanket-ignore `.msh`, `.su2`, `.csv`, `.pdf` or other potentially valuable scientific source formats.

Before a significant simulation, decide how its essential artifacts will be retained. Git-ignored output alone is not archival storage. Keep small publishable manifests, summaries and checksums in tracked case documentation; place large retained artifacts in an appropriate storage location with a retrievable identifier. Choose LFS, releases or scientific storage only when size and retention needs justify it.

For each significant run, record at least:

- A unique run identifier and date; purpose and evidence status.
- Physical conditions, units, reference quantities, model and boundary conditions.
- Exact solver configuration, numerical methods and convergence criteria.
- Mesh source, generation/adaptation settings, element counts and quality indicators; mesh checksum.
- Solver version and Docker image identity/digest, executable/build information where available; Python environment and dependency versions when used.
- Git commit, dirty-state status and any preserved patch needed to identify the actual executed code.
- Exact launch command, process/thread counts, hardware/resource context and relevant execution settings.
- Solver logs, convergence histories, failures, quantities of interest, cost measurements and post-processing commands.
- Retained artifact locations/checksums, data provenance, limitations and verification/validation conclusions.

Never overwrite original observations with transformed data. Record extraction, normalization, filtering and plotting steps. Define and label quantities, units and scientific conventions consistently. Preserve meaningful failures and negative findings.

The unexecuted orchestrator writes `run_manifest.json`, snapshots, `adaptation_history.csv` and per-solve logs/configurations. Source/input/dependency hashes, Git commit/status, UTC time, host Python/NumPy and inspected container image identity are recorded. The baseline snapshot is in `baseline/`; each subiteration NPZ describes its previous mesh/solution, while the new mesh and solve outputs are stored alongside it. Expected solver version is labelled expected; no exact build provenance is invented. Full study-level resource context, retention identifiers, physical conventions, numerical uncertainty and validated convergence evidence still require case documentation and future qualification. This initial manifest is not a complete archival solution.
