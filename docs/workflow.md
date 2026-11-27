# Working and publication workflow

## Repository organization

The seven numbered stage directories organize scientific responsibilities. The baseline-mesh stage now contains real native geometry, a task-specific verifier and one reviewed reference snapshot; other stages remain documentation only. Add implementation, cases and tests only when real inputs and tasks require them.

`AGENTS.md` owns permanent collaboration rules; `scientific-framework.md` owns the scientific scope and methodology; `environment.md` owns tooling observations and environment requirements; `project-status.md` owns current progress, decisions and next actions. The README is the public entry point, not a duplicate project log.

### Stage interactions and shared software

`00_resources` supports all stages. `01_baseline_mesh` owns reusable geometry and baseline mesh definitions. `02_mesh_adaptation` owns the adaptation method and its verification evidence. `03_adaptive_workflow` coordinates those tools and SU2. `04_euler_simulations` and `05_rans_simulations` supply case-specific inputs to the same workflow. `06_experimental_validation` compares approved results with traceable experimental evidence. RANS mesh construction may reuse baseline geometry and the Euler-adapted exterior, subject to adequacy checks.

Numbered directories are scientific stages, not Python package names. When executable code arrives, place reusable adaptation, orchestration and quantitative analysis in one ordinary importable package under root `src/`, with actual tests under root `tests/`; neither is created yet. Scientific-stage documents and cases will refer to those modules. Establish packaging and approved dependencies from the actual imported code, without `sys.path` patches or copies across stages.

Keep dependencies directed: numerical mesh/adaptation code must not import Euler/RANS cases or orchestration; orchestration calls reusable components; cases provide data/configuration; analysis consumes outputs. Share geometry by referencing its authoritative `.geo` definition. ParaView scripts will run through its own application interpreter, while host numerical code uses the designated project venv. No pipeline or circular dependency is implemented.

## Git and GitHub

Use `main` as the primary branch. For now, small verified initialization/documentation changes can go directly to it. Use focused branches and pull requests when changes merit scientific review, involve collaborators or would disrupt the main workflow. Use tags for meaningful reproducible milestones, not every edit.

Before work, inspect status, branch, recent history and remotes. Before synchronization, fetch remote state and inspect divergence; preserve remote contributions. If concurrent changes exist, reconcile non-destructively; request arbitration for substantive scientific conflicts. Never force-push or rewrite history without approval.

Before committing:

1. Inspect every changed file, including untracked content, and the staged diff.
   Before any staging, verify `.PROJECT_CFDvenv/`, `.vscode/` and `local/` are ignored with `git check-ignore`; ensure none are already tracked.
2. Run appropriate checks and record meaningful outcomes. For documentation-only work, check links, consistency, whitespace and publication scope; no invented test suite is required.
3. Stage explicit paths. Review for sensitive information, redistribution rights and unintended large files. Ignore rules are not a secret scanner.
4. Commit one coherent objective with an English Conventional Commit message (`docs:`, `chore:`, `feat:`, `fix:`, `test:` as appropriate).
5. Push ordinary publishable commits to the configured upstream. Verify local/remote commit equality and final status. If access fails, retain local commits and report the blockage accurately.

The repository is intended to be public. Ordinary commits and pushes are authorized; a license grant, destructive operations and uncertain-rights publication still require user approval. Authentication remains in the existing credential mechanism, never in versioned files. Use a repository-local public GitHub noreply address for commit privacy; do not change global identity settings.

## Existing work and reference material

Before importing external code, course documents or datasets, establish origin, authorship, permitted use and redistribution rights. When uncertain, keep them outside the Git working tree and ask the user. Preserve originals and history when available; record justified modifications separately. Link to protected scientific publications instead of redistributing them. Record bibliographic metadata and persistent identifiers only after verification.

Version public reference identifiers, verified metadata and original summaries in `00_resources/`. For protected originals outside the checkout, use an ignored `local/resource-paths.tsv` mapping from a real reference identifier to an external file path. Create entries only when documents actually arrive. The map is local configuration, not part of the public bibliography or an archival copy of the documents.

## Data and generated results

Version source code, native geometry programs, solver configurations, justified parameters, scripts, tests and publishable documentation. Small, important reference meshes or datasets may be tracked after explicit size, provenance and rights review.

Root `runs/` and `outputs/` are ignored locations for generated computations and derived figures/tables, respectively. Future runs should use stage/case/run/cycle identifiers below these roots; exact manifest structure remains to be designed around a real case. Keep inputs and reviewed reference data in their responsible stage, distinguish originals from transformed data, and retain small publishable summaries/manifests there. Do not write generated CFD results into numbered source directories.

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

No run-manifest schema or automation is implemented yet: it will be designed around the first real case and confirmed solver outputs.
