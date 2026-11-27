# Anisotropic Mesh Adaptation for Transonic Airfoil CFD

Scientific project developing and investigating anisotropic triangular mesh adaptation for capturing shocks in transonic airfoil flows, followed by an assessment of viscous effects and boundary-layer resolution.

**Status.** The [NACA 0012 baseline mesh](01_baseline_mesh/README.md) has passed geometric verification. [Mission 02](02_mesh_adaptation/README.md) implements one theory-based adaptation and verifies it on manufactured density fields, including on the baseline mesh. CFD solution verification and experimental comparison are downstream stages.

## Scientific approach

The sequence covers theoretical analysis, numerical implementation and verification, Euler solution/adaptation cycles, transition to RANS with a hybrid boundary-layer mesh, and comparison with experimental literature. An optional study may compare Euler and different RANS mesh strategies using quantitative shock and aerodynamic indicators.

The algorithm currently supports triangles only. Extending it to quadrilaterals or hybrid meshes is outside the main scope. Comparison with the unavailable instructor-provided SU2 adaptation procedure has been abandoned. Independent verification will use controlled tests and convergence studies.

## Project navigation

| File | Responsibility |
| --- | --- |
| [AGENTS.md](AGENTS.md) | Permanent collaboration rules and session entry procedure |
| [Scientific framework](docs/scientific-framework.md) | Scope, roadmap, verification principles and unresolved scientific choices |
| [Workflow](docs/workflow.md) | Git, publication, data provenance and reproducibility practices |
| [Environment](docs/environment.md) | Observed tooling and outstanding environment information |
| [Project status](docs/project-status.md) | Current state, organization decisions and next actions |

## Scientific stages

| Directory | Responsibility |
| --- | --- |
| [00_resources](00_resources/README.md) | References, original notes and publication rights |
| [01_baseline_mesh](01_baseline_mesh/README.md) | Native parametric Gmsh geometry and initial meshes |
| [02_mesh_adaptation](02_mesh_adaptation/README.md) | Triangular adaptation method and verification evidence |
| [03_adaptive_workflow](03_adaptive_workflow/README.md) | Shared automated solve/adapt workflow |
| [04_euler_simulations](04_euler_simulations/README.md) | Euler cases, convergence and shock-capture studies |
| [05_rans_simulations](05_rans_simulations/README.md) | Viscous cases, hybrid boundary-layer meshes and analysis |
| [06_experimental_validation](06_experimental_validation/README.md) | Experimental comparison, provenance and uncertainty |

The baseline stage contains geometry, verification tools and a selected reference snapshot. The adaptation stage owns the theory and verification evidence. Reusable Python code lives in `src/cfd_adapt/` and its tests in `tests/`. Generated runs live outside the versioned source directories. See the [workflow](docs/workflow.md).

## Planned tools

The qualified runtime baseline is host Python 3.14.7, Gmsh `4.15.2-git`, ParaView 6.1.1 and SU2 8.5.0 in the existing Docker container. Host, SU2-wrapper and ParaView Python environments are separate. The AMD64 SU2 image runs under emulation on an ARM64 host; the high-level SU2 Python package lacks NumPy internally, while the low-level wrapper import passes. See [environment qualification](docs/environment.md) for evidence, launch diagnostics and limitations. There is no implemented CFD-case execution procedure yet.

## Reuse and publication

No license has been applied. A public repository alone does not establish permission to reuse its contents. A proposed license for original code is recorded in [project status](docs/project-status.md), pending owner approval and clarification of rights to pre-existing work. Third-party materials require separate provenance and redistribution review; protected publications and teaching materials must not be uploaded without permission.
