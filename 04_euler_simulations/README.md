# Euler simulations

Own approved Euler case inputs, physical/numerical choices, independent verification studies and quantitative shock/accuracy/cost analysis. No Euler case or CFD result was created or executed during Mission 03.

The shared [adaptive workflow](../03_adaptive_workflow/README.md) routes `SOLVER=EULER` here. An omitted solver defaults to Euler and is made explicit in local solve configurations. Each case normally contains:

```text
case_name/
├── baseline.su2
├── config.cfg
├── adaptation.cfg
└── results/
    ├── run_001/
    └── run_002/
```

Case names and inputs are user-owned; no example is generated in the implementation mission. External configurations create uniquely named cases only at future runtime. Existing cases are reused without replacing their inputs, and every run snapshots the actually selected files/dependencies.

Each run contains `baseline/`, `cycle_NN/subiter_NN/`, exact reference/adaptation configurations, a provenance manifest and `adaptation_history.csv`. Each solve retains its mesh, `mesh.png`, local `config.cfg`, `su2.log`, ASCII restart, VTU and CFD history. `adaptation.pvd` animates all completed solutions including the baseline; cycle-specific `cycle_NN.pvd` collections contain that cycle's solutions. PVD indices denote adaptation steps, not physical time. See the main guide for the full tree and parameter format.

Per-case `results/` is ignored; input meshes and configurations remain eligible for deliberate rights/size review and versioning. Root `runs/baseline/` remains reserved for existing baseline generation artifacts, not adaptive Euler runs. Retain publishable study summaries/manifests with case documentation and decide long-term storage before significant computation.

Reference the authoritative geometry and four shared programs without duplication. Distinguish solver-reported iterative convergence, remeshing status, adaptation-cycle stabilization and spatial discretization convergence. No physical conditions or numerical-study criteria are selected by the orchestration example. Independent verification may use controlled cases, conventional refinement and non-adapted baselines; the unavailable instructor procedure remains excluded.
