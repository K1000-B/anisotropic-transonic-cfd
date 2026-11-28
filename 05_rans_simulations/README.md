# RANS simulations

Own viscous/turbulent case definitions, model/wall-treatment choices, boundary-layer and hybrid-mesh specifications, wall-resolution assessment and quantitative scientific analysis. No RANS calculation, hybrid mesh or boundary-layer adaptation was validated or executed during Mission 03.

The shared [adaptive workflow](../03_adaptive_workflow/README.md) routes `SOLVER=RANS` here, using the same case/run convention as Euler:

```text
case_name/
├── baseline.su2
├── config.cfg
├── adaptation.cfg
└── results/
    ├── run_001/
    └── run_002/
```

At future execution, each new run preserves reference/adaptation configurations, dependencies and a manifest, then stores a baseline solve and `cycle_NN/subiter_NN/` results. Each solve retains its mesh, plot, local CFD configuration, solver log, ASCII restart, VTU and convergence history. The global `adaptation.pvd`, cycle-specific `cycle_NN.pvd` and `adaptation_history.csv` record completed steps. Numerical results are ignored under per-case `results/`; case inputs are not automatically excluded from publication review. No example case was created.

The current remesher supports **two-dimensional triangular meshes only**. A RANS destination does not establish compatibility with hybrid meshes containing quadrilateral boundary-layer elements. Such meshes are rejected, never silently triangulated. Hybrid RANS boundary-layer mesh adaptation is not implemented in this mission and belongs to a later project phase.

The longer-term strategy is a justified turbulence model/wall treatment and, initially, an Euler-adapted triangular exterior adjoining a quadrilateral boundary-layer region. Native Gmsh should construct that region where feasible. Verify geometric conformity, transition quality, displaced shocks, viscous gradients and exterior adequacy; reuse is provisional. A purely triangular RANS comparison remains optional and scientifically undecided.

Future evidence must include achieved wall resolution, iterative convergence, targeted mesh sensitivity, aerodynamic quantities and shock–boundary-layer analysis. Reuse authoritative geometry, orchestration and post-processing; separate model effects, numerical uncertainty and experimental validation. The presence of routing code makes no numerical or physical adequacy claim.
