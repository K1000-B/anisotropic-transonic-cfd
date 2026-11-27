# Project status and decisions

Last updated: 2026-11-27.

## Current state

Initialization, runtime qualification and stage organization are complete. The user selected NACA 0012 and authorized baseline mesh generation. The baseline stage now contains native reusable geometry, a parametric triangular C-layout, a specialized verifier and a reviewed reference snapshot (7154 nodes, 13702 triangles). Gmsh coherence, quantitative geometry/topology/quality checks, native SU2 export comparison and visual inspection passed. Default MSH regeneration was byte-identical. No CFD run or physical validation was performed. Details, geometry convention, metrics and limitations are owned by `01_baseline_mesh/README.md`.

Mission 02 implements and verifies one solver-independent adaptation from nodal density using the published project theory. The reference cases include two rotated manufactured transitions and the Mission 01 C-grid. The implementation exports SU2 meshes and ParaView diagnostics. [Mission 02 evidence](../02_mesh_adaptation/verification.md) records geometry checks, residual noncompliance, runtime and achieved node counts. [Mission 03](../03_adaptive_workflow/README.md) owns the subsequent solve/adapt workflow.

Public repository: [K1000-B/anisotropic-transonic-cfd](https://github.com/K1000-B/anisotropic-transonic-cfd). `main` tracks HTTPS `origin/main`. Git history records actual commits and publication history; check local/remote state at each session rather than treating this document as a live synchronization monitor.

## Organization decisions

| Decision | Rationale |
| --- | --- |
| Root instructions plus focused `docs/` files | Durable session entry; separate ownership of rules, science, environment and progress |
| Seven numbered scientific-stage directories with local READMEs | Make the requested scientific sequence explicit without fake programs or data |
| Shared importable code under `src/` and tests under `tests/` | Keep stage names independent from package imports; prevent duplicate algorithms/workflows and circular dependencies |
| Existing designated host venv; no new dependency installation | Respect the supplied environment; record its inventory without inventing project requirements |
| Non-interactive Docker exec with explicit wrapper activation | Avoid terminal interaction and distinguish executable availability from Python initialization |
| Ignored local path settings and external-document mapping | Keep private host paths and protected originals out of public content |
| Root generated run/output locations; reviewed references remain trackable | Separate scientific inputs from large reproducible artifacts without ignoring entire data formats |
| `main`, atomic Conventional Commits and ordinary pushes | Simple workflow; branches/review when useful; no history rewriting |
| No custom skill, CI, automation or license grant | No present need for extra infrastructure; ownership and licensing still require review |

The proposed license for original code remains [MIT](https://choosealicense.com/licenses/mit/), pending explicit approval and ownership/coauthor checks. No license has been applied. Documentation/data or third-party code may require separate terms.

## Completed work and evidence

- Initial foundations were committed and published as `f6f651b` and `5b76b97`; public visibility was verified, including without authentication.
- Inspected the user-designated host Python; checked interpreter identity, package inventory, dependency metadata and imports. No package was added. Details and limitations are in `environment.md` and its linked inventory.
- Inspected Docker image/digest, architecture, mounts, resources and permissions. Started the existing container after explicit user confirmation resolved an automatic approval rejection. The container remains running; its image, configuration and installed software were not changed.
- Verified SU2 8.5.0 help, non-interactive activation, low-level wrapper import, two-rank MPI collectives and shared project read/write access. Diagnostic markers were removed. No solver driver or CFD case was run.
- Qualified Gmsh CLI and ParaView application executables. Both bundled Python tools imported `paraview.simple`; rendering and scientific I/O were not tested.
- Preserved existing editor settings and host venv through explicit ignore rules. Added local machine notes/tool paths outside publication, and numbered-stage READMEs defining responsibilities and inputs/outputs.
- Checked all 15 public candidate files for local links, final newlines, whitespace, unexpectedly large content and selected private-path/credential patterns; no issue was found. Verified ignore rules and absence of already-tracked environment/local files. These bounded checks are not a comprehensive security audit.

## Open limitations and decisions

- AMD64 SU2 runs under emulation on ARM64. Docker has 4 CPUs and about 4.80 GiB RAM; performance/scaling is unmeasured. Resource changes need approval and evidence.
- Low-level `pysu2` imports, but high-level `SU2` import fails due to absent internal NumPy. Use direct executable orchestration initially; request a container change only if a real use requires the high-level tools.
- Runtime checks do not establish mesh/restart format compatibility, a working CFD case, exact upstream source/build provenance, or numerical correctness. Host Python 3.14 compatibility of future dependencies remains to be checked.
- The stage 02 implementation follows the project LaTeX and does not import the previously mentioned external algorithm. Its provenance and rights remain outstanding if it is later supplied. License approval remains outstanding.
- NACA 0012 is selected; its original finite trailing-edge geometry is implemented. Experimental specimen compatibility, domain independence, flow conditions, CFD numerical settings and compute budget remain unresolved. Euler exterior-mesh reuse for RANS remains a strategy to verify.

## Next actions

1. Establish the Mission 03 handoff with a selected CFD case, a real SU2 nodal export and a mesh import check in the qualified container.
2. Set geometry constraints, restart or transfer policy, shared sampling points and separate CFD and remeshing criteria before looping.
3. Assess the baseline geometry convention against candidate transonic experiments and resolve code licensing separately.
