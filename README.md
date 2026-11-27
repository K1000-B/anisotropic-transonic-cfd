# Anisotropic Mesh Adaptation for Transonic Airfoil CFD

Scientific project investigating an existing anisotropic triangular mesh adaptation algorithm for capturing shocks in transonic airfoil flows, followed by an assessment of viscous effects and boundary-layer resolution.

**Status:** project foundations only. No algorithm, geometry, solver configuration, mesh, experimental dataset, or simulation result has been imported or produced. The existing adaptation implementation remains outside this repository pending technical and publication-rights review.

## Scientific approach

The planned sequence is theoretical analysis, implementation review and verification, Euler solution/adaptation cycles, transition to RANS with a hybrid boundary-layer mesh, and comparison with experimental literature. An optional study may compare Euler and different RANS mesh strategies using quantitative shock and aerodynamic indicators.

The algorithm currently supports triangles only. Extending it to quadrilaterals or hybrid meshes is outside the main scope. Comparison with the unavailable instructor-provided SU2 adaptation procedure has been abandoned. Independent verification will use controlled tests and convergence studies.

## Project navigation

| File | Responsibility |
| --- | --- |
| [AGENTS.md](AGENTS.md) | Permanent collaboration rules and session entry procedure |
| [Scientific framework](docs/scientific-framework.md) | Scope, roadmap, verification principles and unresolved scientific choices |
| [Workflow](docs/workflow.md) | Git, publication, data provenance and reproducibility practices |
| [Environment](docs/environment.md) | Observed tooling and outstanding environment information |
| [Project status](docs/project-status.md) | Current state, organization decisions and next actions |

## Planned tools

Python for scientific algorithms, automation and quantitative analysis; native Gmsh `.geo` programs for geometry and initial meshing; SU2 in a user-specified Docker environment for CFD; ParaView for mesh and flow inspection. Versions and dependencies for scientific work are not yet established. There is currently no installation or execution procedure for a CFD case.

## Reuse and publication

No license has been applied. A public repository alone does not establish permission to reuse its contents. A proposed license for original code is recorded in [project status](docs/project-status.md), pending owner approval and clarification of rights to pre-existing work. Third-party materials require separate provenance and redistribution review; protected publications and teaching materials must not be uploaded without permission.
