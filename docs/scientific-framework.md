# Scientific framework

## Objective and scope

Evaluate an existing anisotropic triangular mesh adaptation algorithm for transonic airfoil CFD, especially efficient capture of strong gradients and shocks. Investigate how viscous effects and boundary-layer discretization affect shock predictions, accuracy and computational cost.

The existing implementation is present in `02_mesh_adaptation/`. LaTeX theory was added in commit `e30df92` and removed from the working tree by `1fa3b8b`; it remains in Git history and is not restored or changed here. Numerical correctness remains to be verified; source inspection establishes the operations described below, not successful runtime behavior. Mission 03 adds orchestration without modifying the theoretical optimization problem.

Quadrilateral or hybrid support in the adaptation algorithm is outside the main scope. The unavailable instructor-provided SU2 adaptation procedure will not be used as a comparator.

## Existing implementation and orchestration

The scalar sensor is generic; Mach is the current configuration example, not a selected/validated scientific optimum. Extraction aligns nodal values using restart point IDs and, in the workflow, requires coordinate agreement with the current mesh. No missing physical quantity or transferred solution is reconstructed.

The existing differential reconstruction applies a P1 element-gradient solve and area-weighted nodal recovery, then differentiates and recovers the two gradient components again. It symmetrises the recovered Hessian. Its curvature tensor replaces Hessian eigenvalues by their absolute values bounded below by positive epsilon. Epsilon has scalar-Hessian units and is not a chord-relative length.

The analytical Lp formulation constructs the optimal metric using lumped nodal areas and the prescribed target complexity. Spectral bounds impose minimum/maximum metric sizes and maximum principal size ratio; a scalar scaling with bounding restores complexity. The source implementation owns these formulas. The workflow supplies `HMINS[c] * CHORD` and `HMAX * CHORD`, integer-valued positive p, unchanged epsilon and dimensionless edge thresholds. It reuses the existing complexity-feasibility logic before CFD.

`BUDGET[c]` is an indicative vertex budget converted to `C_star = sqrt(3)/2 * BUDGET[c]`. The final vertex count is not imposed exactly, and the target is not adjusted from its discrepancy. The prescribed metric's restored complexity does not establish the remesher's geometric conformity or convergence.

During one remeshing operation, the continuous P1 metric on the source mesh remains fixed. The existing splits, collapses, swaps and relocations use that field; `REMESH_MAX_ITER` limits their internal passes. Every workflow subiteration extracts the latest solved scalar and recomputes the metric, followed by a new CFD solve on the adapted mesh. `SUBITER` counts these complete steps within a cycle; cycle parameters remain fixed while mesh and solution evolve.

Source-domain robustness safeguards use triangle convexity to integrate an edge contained in one source triangle directly, with the same barycentric tolerance/roundoff cleanup as metric point evaluation. Other edges retain the existing triangle-interval integration. Length and midpoint evaluations share a canonical segment ordering, and every changed candidate edge is checked for source-metric coverage before a local edit is committed. These are implementation safeguards, not a changed Hessian/metric optimization, relaxed domain constraint or extrapolation rule. Their numerical behavior remains to be verified separately.

Airfoil/farfield vertices and non-outlet frozen edges remain protected by the existing policy. Outlet-only vertices can adapt on the existing straight outlet geometry, with markers retained. The workflow adds no boundary reconstruction, external remesher, exact-node-count loop or solution transfer. A valid stationary/iteration-limited mesh may proceed with an explicit diagnostic; this is not a claim of metric convergence.

All solves, including baseline and the final adapted mesh, use fresh reference-flow initialisation. Physical/numerical CFD choices are preserved from the user reference configuration. Solver-reported iterative convergence, remeshing convergence, cycle stabilization, spatial convergence and model validation remain distinct. No simulation or numerical test was executed in Mission 03.

## Roadmap

1. **Theoretical framework.** Explain the mathematical and numerical foundations actually used by the existing algorithm, including assumptions and limitations.
2. **Numerical implementation.** The elementary implementation and orchestration are present; inspect/document them and verify in a separate authorized testing phase. Develop meaningful tests and make improvements only when justified; preserve original work and provenance.
3. **Euler verification study.** Select an airfoil and primarily one Mach number consistent with available experimental evidence. Justify physical conditions, geometry, initial mesh and numerical settings. Assess iterative convergence, successive solve/adapt cycles and mesh convergence. Quantify shock capture, accuracy and cost against appropriate non-adapted or conventionally refined meshes and controlled tests.
4. **Euler-to-RANS transition.** Initially construct a quadrilateral boundary-layer region and use the Euler-adapted triangular mesh as a starting point for the exterior. Verify geometric conformity, element quality and the transition. Justify the turbulence model, wall treatment, target and achieved wall resolution, numerical settings and iterative convergence. Establish adequacy for viscous effects and shock–boundary-layer interaction. Propose targeted mesh-sensitivity checks when needed rather than automatically repeating the entire Euler study.
5. **Experimental comparison.** Identify primary references and accessible datasets. Check comparability of physical conditions and definitions; compare RANS predictions quantitatively and discuss experimental, numerical and modelling uncertainties. Record data and redistribution provenance.
6. **Optional investigation.** Compare Euler on adapted triangles, RANS on a hybrid mesh, and potentially RANS on adapted triangles. A less restrictive near-wall resolution is admissible only if compatible with the actual model and wall treatment, with limitations documented. Separate physical-model effects from discretization effects as far as possible.
7. **Conclusions.** Assess adaptation efficiency and criteria, shock accuracy, viscous and boundary-layer effects, accuracy/cost tradeoffs, limitations and justified improvements.

## Evidence and verification principles

- **Code verification:** test mathematical operations and implementation against exact or controlled cases where available. Establish tests from the actual algorithm after import.
- **Solution verification:** estimate discretization and iterative errors for the computed quantities of interest. Separate iterative convergence from mesh convergence.
- **Adaptation convergence:** stabilization of mesh/cycle indicators is a separate stopping concept. It establishes neither sufficient convergence of each CFD solve nor spatial discretization convergence. Define and assess all three independently.
- **Model validation:** compare physical predictions with suitable experimental evidence, accounting for differences and uncertainty. Euler results are not a viscous experimental validation by default.
- Quantify shock location, strength and numerical thickness with explicit extraction definitions; pressure distributions, lift, drag, boundary-layer characteristics, shock–boundary-layer interaction, mesh complexity and computational cost may also be relevant.
- Specify dimensional or nondimensional conventions and normalization before comparisons. Distinguish pressure and viscous drag contributions where relevant.
- Compare adaptation strategies at meaningful accuracy/cost levels; include adaptation overhead in cost accounting. A finer-looking mesh or sharper-looking plot does not establish improved accuracy.
- Do not assume smooth-solution error orders or uncritically apply extrapolation-based uncertainty estimates to shock-dominated, anisotropically adapted mesh sequences. Assess whether the chosen convergence method is applicable.
- An Euler-adapted mesh may miss shifted shocks, wall gradients or other viscous features. Hybrid conformity and exterior adequacy must both be demonstrated.
- Track residuals together with quantities of interest; stopping criteria and acceptable sensitivities remain to be defined and justified.

## Unresolved scientific inputs

The user selected NACA 0012 on 2026-11-27. The baseline uses its original finite-thickness trailing edge; the meshing stage documents the exact definition and initial domain/sizing choices. These do not establish domain independence or agreement with a future experimental specimen. No Mach number, incidence, Reynolds number, thermodynamic state, physical reference quantities, turbulence model, wall treatment, wall-resolution target, flux scheme, limiter, adaptation variable, stopping tolerance or computational budget has been selected. No reference experimental dataset has been approved. These choices require evidence and, for important scientific decisions, user arbitration.
