# Anisotropic Mesh Adaptation for Transonic Airfoil CFD

Scientific project investigating Hessian-based anisotropic triangular mesh adaptation for transonic airfoil shock capture, followed by an assessment of viscous effects and boundary-layer resolution.

**Status:** a native parametric NACA 0012 baseline mesh was generated and geometrically verified (7154 nodes, 13702 triangles). Three elementary Python programs and an automated solve/extract/adapt orchestrator are present. Mission 03 implementation and static source review are complete; the workflow has not been executed. No CFD result, numerical validation of the complete workflow or experimental validation is established. See [project status](docs/project-status.md) for the evidence boundaries.

## Scientific approach

The sequence is theoretical analysis and code verification, Euler solution/adaptation and spatial-sensitivity studies, transition to RANS with justified wall resolution, and quantitative comparison with experimental literature. The existing adaptation method supports 2D triangles only. Hybrid/quadrilateral boundary-layer mesh adaptation is a later phase; a RANS output route does not establish support for it. Reuse of an Euler-adapted exterior remains a strategy requiring verification. Comparison with the unavailable instructor SU2 procedure is excluded; verification is independent.

## Repository structure

| Directory | Responsibility and current contents |
| --- | --- |
| [00_resources](00_resources/README.md) | Reference metadata, original notes and publication-rights tracking |
| [01_baseline_mesh](01_baseline_mesh/README.md) | Reusable native Gmsh geometry, specialized verifier and reviewed baseline snapshot |
| [02_mesh_adaptation](02_mesh_adaptation/README.md) | Authoritative extraction, triangular adaptation and Docker CFD runner programs; method documentation |
| [03_adaptive_workflow](03_adaptive_workflow/README.md) | One procedural `adaptive_workflow.py` and its user guide |
| [04_euler_simulations](04_euler_simulations/README.md) | Euler case inputs and ignored per-case `results/run_NNN/` |
| [05_rans_simulations](05_rans_simulations/README.md) | RANS case inputs and the equivalent ignored results organization |
| [06_experimental_validation](06_experimental_validation/README.md) | Experimental comparison, provenance and uncertainty; no dataset/result yet |
| `docs/` | Durable scientific, environment, workflow and status documentation |
| `local/` | Ignored machine settings and preserved historical local archives |
| `runs/`, `outputs/` | Ignored baseline/other generated artifacts and derived report outputs |

The four-program architecture is `adaptive_workflow.py` calling `run_simulation.py`, `extract_solution.py` and `mesh_adaptation.py` directly. Host Python performs orchestration, extraction, adaptation, plotting and CSV/PVD management. Only SU2 CFD resolution runs inside the existing Docker container. The elementary sources remain in stage 02; they are not duplicated between cases or moved to a new package.

Each simulation case normally contains `baseline.su2`, `config.cfg`, `adaptation.cfg` and its own `results/run_001/`, `run_002/`, etc. The reference CFD configuration identifies Euler or RANS automatically. External/differently located inputs are snapshotted without altering originals. Each run retains baseline and cycle/subiteration meshes, configurations, restart/VTU/history outputs, plots, provenance and global/cycle PVD collections. No example case or numerical output was created in Mission 03.

## Future execution

Use the existing `.PROJECT_CFDvenv` and verify interpreter identity according to [environment.md](docs/environment.md) before executing project Python. In a separate user testing session, from the checkout root:

```bash
source .PROJECT_CFDvenv/bin/activate
python 03_adaptive_workflow/adaptive_workflow.py
```

The program asks for the mesh, reference CFD configuration and adaptation configuration one at a time. Optional non-interactive execution accepts exactly those three positional paths. The complete parameter format, output tree, dependency handling, error behavior and pending runtime checks are in the [workflow user guide](03_adaptive_workflow/README.md). These commands were not executed during Mission 03.

The intended environment is native macOS project Python with NumPy/Matplotlib, the existing `su2-docker` container and SU2 8.5.0. Historical runtime qualification does not establish integration of this new source. AMD64 SU2 executes under emulation on ARM64; high-level container `SU2` Python tools lack internal NumPy, and are not used by this workflow. Nothing was installed or changed in the container.

## Project documentation

| File | Responsibility |
| --- | --- |
| [AGENTS.md](AGENTS.md) | Collaboration, session entry, scope and publication rules |
| [Scientific framework](docs/scientific-framework.md) | Implemented method, scientific roadmap and unresolved choices |
| [Workflow](docs/workflow.md) | Execution architecture, Git, provenance and result retention |
| [Environment](docs/environment.md) | Historical observations, intended execution domains and limitations |
| [Project status](docs/project-status.md) | Current implementation and verification status, next actions |

## Reuse and publication

No license has been applied. A public repository alone does not establish permission to reuse its contents. The proposed license for original code remains subject to owner approval and rights checks. Third-party material requires separate provenance/redistribution review; protected scientific papers, teaching materials and datasets must not be uploaded without permission. Historical theory and local archives were preserved without restoring files removed by earlier commits.
