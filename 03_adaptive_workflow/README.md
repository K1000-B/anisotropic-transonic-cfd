# Automated adaptive CFD workflow

Own the shared solve/adapt cycle contract and orchestration documentation. No routine is implemented yet. Reuse one implementation for all cases rather than duplicating scripts.

Future inputs are geometry/mesh references, reviewed case configuration, adaptation controls and local tool/path settings. The workflow will coordinate Gmsh, host Python and SU2 Docker: prepare inputs, solve, assess iterative convergence, extract fields, adapt, check quality, compare quantities and decide whether to continue.

Outputs will include cycle-indexed meshes, configurations, histories, quantities of interest, costs, error logs and provenance under root `runs/`. Failure handling and restart checkpoints must be designed around real solver behavior. A failed cycle must remain identifiable and must not silently produce an accepted next mesh.

Treat iterative CFD convergence, mesh adaptation convergence and spatial discretization convergence separately. A stable adaptation stopping criterion is not a mesh-error estimate. The verified runtime launch method is in [environment qualification](../docs/environment.md).
