# Project status and decisions

Last updated: 2026-11-28.

## Current state

Initialization, historical runtime qualification and baseline generation are complete. The user selected NACA 0012. Its stage contains native reusable geometry, a parametric triangular C-layout, a specialized verifier and a reviewed reference snapshot (7154 nodes, 13702 triangles). Gmsh coherence, quantitative geometry/topology/quality checks, native SU2 export comparison and visual inspection passed; default MSH regeneration was byte-identical. These previously established results are preserved. Details and limitations are owned by `01_baseline_mesh/README.md`.

Three elementary Python programs are present in stage 02 (`1fa3b8b`). That commit removed the LaTeX files previously added by `e30df92`; their historical state is preserved without restoration. Mission 03 implements `03_adaptive_workflow/adaptive_workflow.py`, minimal diagnostics/runner compatibility changes, per-case results ignore rules and consistent documentation. Review was static only: no Python interpreter/import/compilation, test, fixture, Docker/SU2 command, mesh/plot/result generation or installation was performed during this mission. No CFD solution or physical validation is established.

| Component | Status |
| --- | --- |
| Baseline mesh generation | Previously generated and geometrically verified; evidence retained |
| Existing extraction implementation | Present; source unchanged in Mission 03 |
| Existing mesh adaptation implementation | Present; optional observation-only diagnostics added |
| Existing Docker CFD runner | Present; minimal launch/log/mount/interruption compatibility changes |
| Automated adaptive workflow | Implemented, pending execution |
| Repository documentation | Updated for actual interfaces and architecture |
| Static source code review | Completed for Mission 03 changes; no execution-based checks |
| Workflow runtime testing | Not performed in Mission 03 |
| Docker integration testing | Not performed in Mission 03; historical diagnostics below remain valid within their original scope |
| CFD integration testing | Not performed |
| Numerical validation of the complete workflow | Pending |
| Performance assessment | Pending |
| Hybrid RANS mesh adaptation | Not implemented in this mission |
| Experimental validation | No approved dataset or validation result |

Public repository: [K1000-B/anisotropic-transonic-cfd](https://github.com/K1000-B/anisotropic-transonic-cfd). `main` remains the primary branch; Mission 03 continues the existing `codex/02-mesh-adaptation` branch and its HTTPS origin upstream. Git history records commits/publication; check current synchronization each session rather than treating this document as a live monitor.

## Organization decisions

| Decision | Rationale |
| --- | --- |
| Root instructions plus focused `docs/` files | Durable session entry; separate ownership of rules, science, environment and progress |
| Seven numbered scientific-stage directories with local READMEs | Make the requested scientific sequence explicit without fake programs or data |
| Preserve the three elementary stage-02 programs and one procedural stage-03 orchestrator | Mission 03 explicitly supersedes the earlier proposed `src/` relocation; direct function imports without duplicated algorithms/classes |
| Existing designated host venv; no new dependency installation | Respect the supplied environment; record its inventory without inventing project requirements |
| Non-interactive Docker exec with explicit wrapper activation | Avoid terminal interaction and distinguish executable availability from Python initialization |
| Ignored local path settings and external-document mapping | Keep private host paths and protected originals out of public content |
| Per-case adaptive `results/run_NNN/`; root baseline runs/report outputs retained | Separate Euler/RANS calculations and source inputs; reviewed references remain trackable without blanket format ignores |
| `main`, atomic Conventional Commits and ordinary pushes | Simple workflow; branches/review when useful; no history rewriting |
| No custom skill, CI, automation or license grant | No present need for extra infrastructure; ownership and licensing still require review |

The proposed license for original code remains [MIT](https://choosealicense.com/licenses/mit/), pending explicit approval and ownership/coauthor checks. No license has been applied. Documentation/data or third-party code may require separate terms.

## Completed work and evidence

- Initial foundations were committed and published as `f6f651b` and `5b76b97`; public visibility was verified, including without authentication.
- Inspected the user-designated host Python; checked interpreter identity, package inventory, dependency metadata and imports. No package was added. Details and limitations are in `environment.md` and its linked inventory.
- Inspected Docker image/digest, architecture, mounts, resources and permissions. Started the existing container after explicit user confirmation resolved an automatic approval rejection. It was left running at that assessment; Mission 03 did not check its current state. Its image, configuration and installed software were not changed.
- Verified SU2 8.5.0 help, non-interactive activation, low-level wrapper import, two-rank MPI collectives and shared project read/write access. Diagnostic markers were removed. No solver driver or CFD case was run.
- Qualified Gmsh CLI and ParaView application executables. Both bundled Python tools imported `paraview.simple`; rendering and scientific I/O were not tested.
- Preserved existing editor settings and host venv through explicit ignore rules. Added local machine notes/tool paths outside publication, and numbered-stage READMEs defining responsibilities and inputs/outputs.
- Checked all 15 public candidate files for local links, final newlines, whitespace, unexpectedly large content and selected private-path/credential patterns; no issue was found. Verified ignore rules and absence of already-tracked environment/local files. These bounded checks are not a comprehensive security audit.
- Mission 03 source review covers direct imports/signatures, parameter parsing/propagation, sequential prompts, solver/case routing, fresh initialisation, latest mesh/solution propagation, final solve, protected-boundary compatibility, plots, CSV/PVD, snapshots, errors and interruption handling. Upstream SU2 8.5.0 template/output/driver source was inspected for execution/output conventions. Installed-build behavior remains pending; the [workflow guide](../03_adaptive_workflow/README.md) owns the detailed contracts and limitations.

## Open limitations and decisions

- AMD64 SU2 runs under emulation on ARM64. Docker has 4 CPUs and about 4.80 GiB RAM; performance/scaling is unmeasured. Resource changes need approval and evidence.
- Low-level `pysu2` imports, but high-level `SU2` import fails due to absent internal NumPy. Use direct executable orchestration initially; request a container change only if a real use requires the high-level tools.
- Runtime checks do not establish mesh/restart format compatibility, a working CFD case, exact upstream source/build provenance, or numerical correctness. Host Python 3.14 compatibility of future dependencies remains to be checked.
- Existing algorithm and historical theory are already in Git; their full provenance/redistribution and licensing review remains outstanding. Mission 03 publishes new orchestration/documentation and minimal compatibility edits, not new external scientific resources. License approval remains outstanding; no license was applied.
- NACA 0012 is selected; its original finite trailing-edge geometry is implemented. Experimental specimen compatibility, domain independence, flow conditions, CFD numerical settings and compute budget remain unresolved. Euler exterior-mesh reuse for RANS remains a strategy to verify.

## Next actions

1. Conduct a separate testing and debugging session with the user: first review a real case, dependencies, physical assumptions and compute budget, then qualify Python/imports, configuration/output naming and Docker integration within that session's authorization. Do not execute those checks as part of Mission 03.
2. Verify the elementary numerical operations independently and assess complete workflow behavior/performance before drawing solution or adaptation-convergence conclusions. Confirm large-mesh feasibility, generic sensor output groups and interruption cleanup.
3. Complete provenance/rights and license arbitration before applying a license. The high-level container Python package is not needed by the implemented executable runner; no repair is currently authorized.
4. Check baseline geometry against candidate transonic experiments before selecting operating conditions; plan solution verification and uncertainty assessment separately. Hybrid RANS construction/adaptation and experimental validation remain later scientific phases.
