# Project collaboration instructions

## Session entry and continuity

1. Read this file, `docs/project-status.md` and `docs/environment.md` before acting.
2. Read `docs/scientific-framework.md` and `docs/workflow.md` before relevant scientific, data, Git or publication work.
3. Inspect `git status --short --branch`, current branch, recent commits and remotes. Preserve all existing contributions. Check remote state before synchronization when network access is available.
4. Use the documents above as the durable source of truth, not chat history or implicit memory. Update the responsible document when decisions, restrictions or project state change. Do not duplicate logs across files.

## Communication and scientific judgment

- Discuss the project with the user in French. All deliverables must be in English: code, comments, documentation, figures, tables, captions, commit messages, file and directory names.
- Assume advanced engineering knowledge. Work rigorously, critically and reproducibly; challenge unsupported assumptions and justify significant physical, mathematical, numerical and software choices.
- Prefer simple, modular, maintainable solutions. Avoid speculative infrastructure, unnecessary abstractions and duplicate functionality.
- Never invent results, experimental data, references, software options or successful actions. Distinguish facts, hypotheses, estimates, planned work and verified results. Report limitations, failures and inconsistencies explicitly.
- Do not silently change approved assumptions or scientific decisions. Keep the global objectives in view.
- Decide routine technical matters autonomously. Ask the user about substantial scientific or scope changes, important physical assumptions, significant unexpected computation costs, software-environment changes, data risks and uncertain publication rights.
- Work incrementally; do not implement later scientific phases before they are requested.

## Scientific constraints

- The existing adaptation algorithm supports triangles only; preserve and examine it before proposing justified changes. Do not rewrite it wholesale or extend element support by default.
- The comparison with the unavailable instructor SU2 adaptation procedure is definitively excluded. Use independent verification.
- Treat reuse of an Euler-adapted exterior mesh for RANS as an initial strategy requiring verification, not as established adequacy.
- Distinguish code verification, solution verification, iterative CFD convergence, mesh adaptation convergence, spatial discretization convergence, model validation, numerical uncertainty and physical modelling uncertainty. Stable adaptation cycles do not establish spatial convergence.
- Base conclusions on quantitative, reproducible evidence and appropriate quantities of interest. A solver exit status or visual comparison alone is insufficient.
- Preserve original data; document transformations, units, conventions, provenance, uncertainties and relevant negative results. Prefer primary scientific sources and official documentation.

## Tool and environment constraints

- Use `.PROJECT_CFDvenv/bin/python` for host project Python work. Before every execution or installation, verify `sys.executable` and `sys.prefix` identify that environment; do not rely only on an activated prompt. The user authorizes justified dependencies in this environment, but no global or anticipatory installations. Do not create another environment without approval. Record actual requirements when needed; the observed package inventory is not a requirements specification.
- Keep `.PROJECT_CFDvenv/`, `.vscode/`, `local/`, `runs/` and `outputs/` ignored; verify ignore behavior before staging. Never publish machine-specific absolute paths. Use repository-relative paths and ignored local tool/path configuration.
- The reference solver container is `su2-docker`; starting it and running non-CFD qualification diagnostics are authorized. Do not change its image, configuration, mounts or software. Follow the verified non-interactive activation in `docs/environment.md`; an interactive shell and `docker exec` have different Python initialization. Keep host Python, container Python and ParaView's bundled interpreter separate.
- Read the documented limitations before using SU2 Python tools: `pysu2` imports, but the high-level `SU2` package fails because NumPy is absent internally. No container dependency change is authorized. AMD64 execution on the ARM64 Docker host requires emulation; do not assume native timing or parallel scaling.
- Do not modify Docker, installed software, global Codex settings, other projects or user-wide configuration without approval. Prefer repository-local mechanisms.
- Use native, readable, English-commented, parametric Gmsh `.geo` programs for operations its language can express. Avoid Python scripts that merely enumerate geometry points or reproduce Gmsh features. Use Python when it adds real algorithmic or automation value.
- Verify Gmsh syntax against official documentation when needed. Verify every SU2 option, syntax and behavior against the actual solver version; do not copy unexplained tutorial settings.
- Keep SU2 `.cfg` files and launch commands documented and reproducible. Use ParaView state files or scripting when useful, and Python for reproducible quantitative analysis and report-quality figures.
- During foundation and environment qualification work: no CFD runs, generated meshes, scientific algorithms, arbitrary physical/numerical parameters, mock programs or demonstration results.
- Numbered directories organize scientific stages, not Python import namespaces. The three authoritative elementary programs remain in `02_mesh_adaptation/`; the single procedural orchestrator is `03_adaptive_workflow/adaptive_workflow.py`. Preserve their mathematical method and readable direct-call architecture; do not add packages, classes or copies without a real need. Follow `docs/workflow.md` for ownership and protected-resource references.
- Euler/RANS case inputs belong in `04_euler_simulations/` and `05_rans_simulations/`; generated adaptive outputs belong in each case's ignored `results/run_NNN/`. Root `runs/` retains baseline-generation and other documented outputs. Keep configurations, original inputs and generated outputs distinct.
- Mission 03 is implementation-only: source inspection, Markdown updates and authorized Git operations are permitted; no Python execution/import/compilation, tests, fixtures, Docker/SU2 commands, generated meshes/results or dependency changes. This is a mission-specific restriction; a later user testing/debugging session is separate. Preserve `local/mission-02-abandonment-2026-11-28/` and the existing LaTeX theory unchanged.

## Git and public publication

- Git and GitHub are mandatory. Review differences and perform appropriate checks before each atomic commit; commit each coherent verified stage and ordinarily push publishable changes without asking again.
- Use English Conventional Commit messages, for example `docs: clarify verification strategy` or `feat: add mesh reader`.
- Preserve others' edits. Never force-push, rewrite history, perform destructive Git operations or arbitrarily resolve scientific conflicts without authorization.
- Review staged paths and contents for secrets, personal information, local paths, rights and large generated files before publication. Never stage everything indiscriminately.
- Keep uncertain-rights imports outside the repository until reviewed. `local/` is only a convenience ignore rule, not a confidentiality guarantee. Do not upload protected papers, teaching documents or experimental datasets without established redistribution rights.
- Do not apply a license without user approval. Do not introduce LFS, releases, external storage, CI, custom skills or scheduled automations without an actual need.
- Maintain run provenance and source/result separation according to `docs/workflow.md`. Do not indiscriminately ignore all mesh or data formats: selected reviewed reference artifacts may need versioning.

## Task completion

- Verify work with checks appropriate to its impact; do not equate written code with completion.
- Update durable state and decisions where needed, inspect the final diff, commit and synchronize publishable changes. If blocked, continue independent work and report the exact remaining action.
- Give a concise French report distinguishing implemented, tested, verified, physically validated, documented, committed and published work. Only claim actions actually executed.
