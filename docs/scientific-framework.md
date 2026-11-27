# Scientific framework

## Objective and scope

Evaluate an existing anisotropic triangular mesh adaptation algorithm for transonic airfoil CFD, especially efficient capture of strong gradients and shocks. Investigate how viscous effects and boundary-layer discretization affect shock predictions, accuracy and computational cost.

The existing algorithm and supporting work are external to this repository. Its mathematical foundations and actual capabilities must be established from the supplied implementation and documents. Metric-based adaptation, Riemannian metric tensors, Hessian reconstruction, interpolation-error control, complexity optimization, metric intersection and regularization are candidate topics, not confirmed implementation features.

Quadrilateral or hybrid support in the adaptation algorithm is outside the main scope. The unavailable instructor-provided SU2 adaptation procedure will not be used as a comparator.

## Roadmap

1. **Theoretical framework.** Explain the mathematical and numerical foundations actually used by the existing algorithm, including assumptions and limitations.
2. **Numerical implementation.** Inspect, document and verify the existing triangular algorithm. Develop meaningful tests and make improvements only when justified; preserve original work and provenance.
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

No airfoil, Mach number, incidence, Reynolds number, thermodynamic state, domain size, reference quantities, turbulence model, wall treatment, wall-resolution target, flux scheme, limiter, adaptation variable, stopping tolerance or computational budget has been selected. No reference dataset has been approved. These choices require evidence and, for important scientific decisions, user arbitration.
