# RANS simulations

Own viscous/turbulent case definitions, boundary-layer and hybrid-mesh specifications, wall-resolution assessment and scientific analysis. No RANS case, hybrid mesh or result exists yet.

Future inputs are a justified turbulence model and wall treatment, physical conditions, and initially an Euler-adapted triangular exterior mesh. Native Gmsh definitions should construct the quadrilateral boundary-layer region where feasible. Check conformity, transition quality and suitability for displaced shocks and viscous gradients; exterior reuse is provisional.

Outputs will include achieved wall resolution, iterative convergence, targeted mesh sensitivities, aerodynamic quantities and shock–boundary-layer analysis. Reuse shared orchestration and analysis. Store generated data under root `runs/` and `outputs/`. A purely triangular RANS comparison remains optional and undecided.
