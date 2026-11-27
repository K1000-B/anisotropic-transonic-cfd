# Working and publication workflow

## Repository organization

Keep the current repository documentation-only until real inputs arrive. Add source, geometry, cases, tests and data directories when they have actual content and a clear responsibility. Do not introduce placeholder programs, artificial datasets or an unrequested package/build system.

`AGENTS.md` owns permanent collaboration rules; `scientific-framework.md` owns the scientific scope and methodology; `environment.md` owns tooling observations and environment requirements; `project-status.md` owns current progress, decisions and next actions. The README is the public entry point, not a duplicate project log.

## Git and GitHub

Use `main` as the primary branch. For now, small verified initialization/documentation changes can go directly to it. Use focused branches and pull requests when changes merit scientific review, involve collaborators or would disrupt the main workflow. Use tags for meaningful reproducible milestones, not every edit.

Before work, inspect status, branch, recent history and remotes. Before synchronization, fetch remote state and inspect divergence; preserve remote contributions. If concurrent changes exist, reconcile non-destructively; request arbitration for substantive scientific conflicts. Never force-push or rewrite history without approval.

Before committing:

1. Inspect every changed file, including untracked content, and the staged diff.
2. Run appropriate checks and record meaningful outcomes. For documentation-only work, check links, consistency, whitespace and publication scope; no invented test suite is required.
3. Stage explicit paths. Review for sensitive information, redistribution rights and unintended large files. Ignore rules are not a secret scanner.
4. Commit one coherent objective with an English Conventional Commit message (`docs:`, `chore:`, `feat:`, `fix:`, `test:` as appropriate).
5. Push ordinary publishable commits to the configured upstream. Verify local/remote commit equality and final status. If access fails, retain local commits and report the blockage accurately.

The repository is intended to be public. Ordinary commits and pushes are authorized; a license grant, destructive operations and uncertain-rights publication still require user approval. Authentication remains in the existing credential mechanism, never in versioned files. Use a repository-local public GitHub noreply address for commit privacy; do not change global identity settings.

## Existing work and reference material

Before importing external code, course documents or datasets, establish origin, authorship, permitted use and redistribution rights. When uncertain, keep them outside the Git working tree and ask the user. Preserve originals and history when available; record justified modifications separately. Link to protected scientific publications instead of redistributing them. Record bibliographic metadata and persistent identifiers only after verification.

## Data and generated results

Version source code, native geometry programs, solver configurations, justified parameters, scripts, tests and publishable documentation. Small, important reference meshes or datasets may be tracked after explicit size, provenance and rights review.

Root `runs/` and `outputs/` are ignored locations for generated computations. Root `local/` is ignored for local notes or scratch material, but confidential/uncertain-rights originals should remain outside the repository. These directories are not created until needed. Do not blanket-ignore `.msh`, `.su2`, `.csv`, `.pdf` or other potentially valuable scientific source formats.

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
