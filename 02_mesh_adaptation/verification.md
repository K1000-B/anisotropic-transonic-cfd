# Numerical verification

The [reproduction instructions](README.md) run the automated suite and the three manufactured-field cases. The [specification](mathematical-specification.md) maps individual checks to theory labels and implementation functions.

## Cases

- The square cases use a deterministic, slightly perturbed 16-by-16 cell grid and a tanh density transition through the centre. The two normal directions make angles of 19 and 67 degrees to the x-axis.
- The Mission 01 case uses the reference NACA 0012 C-grid with its original finite trailing edge and a tanh transition at 23 degrees through `(1.7, 0.2)` in chord units.
- The quadratic study compares reconstructed and exact Hessians on regular and jittered meshes, with separate interior and boundary statistics.

Verification checks each SU2 export by rereading coordinates, connectivity and markers. The C-grid check retains the polygonal airfoil, outer boundary and embedded wake. The manufactured inputs provide code-verification evidence. Flow-solution verification and physical validation belong to the later CFD stages.

## Component checks

The 25 automated tests cover signed triangle areas, connected vertex fans, simple boundaries, lumped areas, exact constant and affine gradients, double-P1 Hessian recovery, rotated and indefinite spectra, principal directions, size and aspect limits, infeasible complexity, saturation endpoints and bisection. The edge tests include a multi-face piecewise-affine metric with a closed integral, traversal reversal, and traversal along and across a hole boundary. The operation tests cover each local edit, marker inheritance, protected vertices and edges, invalid trial rejection, and distinct termination reasons. I/O tests cover the Mission 01 mesh round trip, reordered solution rows, geometry/ID disagreement, ambiguous matches, missing data, ParaView-format arrays and recovery from an interrupted output write.

The quadratic error uses the Frobenius norm, with nodal areas weighting the domain RMS. The interior sample excludes two mesh spacings near each side. On grids of 12, 24 and 48 cells per side, regular-grid RMS decreased from 0.2447 to 0.1726 to 0.1219, and jittered-grid RMS decreased from 0.2616 to 0.1963 to 0.1529. The regular-grid interior errors were at round-off. On the jittered grids, interior means stayed between 0.078 and 0.086. Boundary means stayed between 0.462 and 0.493. These results show substantial local boundary bias and persistent error under cell-size-proportional jitter, despite a decrease in the domain-wide RMS.

## Performance method

The Cloud run used CPython 3.12.14, NumPy 2.5.5 and Matplotlib 3.11.2 on Linux x86-64, with Intel Xeon Platinum 8573C host CPUs and one requested OpenBLAS thread. Timings use `perf_counter` and cover the adaptation call, excluding input reading, diagnostic exports and plotting. These are single-run observations on this environment.

The four-sweep square-19 development profile found 25,891 calls to the exact edge integrator, mainly from candidate relocation. Batching the barycentric interval endpoints and projecting the nodal tensors onto the edge direction before interpolation reduced profiled integrator time from 9.62 s to 7.84 s and profiled run time from 14.42 s to 12.46 s. The before/after adapted SU2 hashes and all four accepted-operation counts matched. This comparison predates the later subdivision quality guard. The final eight-sweep timings below include the guard.

## Final adaptation cases

Run date 2026-11-27. [Machine-readable evidence](verification-results.json) records the source revision, case inputs, metrics, operation counts, hashes and timings. All three cases use $p=1$, the edge length band $[0.62,1.45]$, and an eight-pass budget. The square tests use $\varepsilon=0.05$, size limits $[0.023,0.35]$ and $A_{\max}=7$. The C-grid uses $\varepsilon=5\times10^{-4}$, size limits $[0.01,4]$ in chord units and $A_{\max}=12$. These are manufactured-case controls.

| Case | Source vertices → adapted vertices (target) | Target discrete complexity | Relative complexity error after bounds | Metric mean-ratio median, source → adapted | Out-of-band edges, source → adapted | Adaptation wall time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Square, 19° | 289 → 409 (280) | 242.48711306 | $3.22\times10^{-12}$ | 0.319 → 0.942 | 511 → 40 | 8.84 s |
| Square, 67° | 289 → 422 (280) | 242.48711306 | $1.99\times10^{-11}$ | 0.292 → 0.939 | 504 → 30 | 10.04 s |
| Mission 01 C-grid, 23° | 7,154 → 4,994 (4,200) | 3637.30669589 | $5.67\times10^{-11}$ | 0.219 → 0.972 | 20,715 → 1,332 | 104.74 s |

All runs terminated at the pass limit. Remaining short/long edge counts were 40/0, 30/0 and 1,230/102, with 38, 30 and 595 classified as constrained violations. Median metric length moved from 1.309 to 0.896, 1.408 to 0.889 and 0.026 to 1.019, respectively. The target-to-actual node discrepancies are +46%, +51% and +19%. The complexity identity uses lumped quadrature on the source mesh, while the conversion from complexity to nodes assumes unit, equilateral metric triangles and a relatively small boundary contribution.

| Case | Accepted splits | Accepted collapses | Accepted swaps | Accepted relocations | Final metric quality minimum | Physical quality median, source → adapted | Final physical minimum |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Square, 19° | 430 | 310 | 748 | 1,118 | 0.5134 | 0.855 → 0.414 | 0.1535 |
| Square, 67° | 504 | 371 | 761 | 1,101 | 0.2944 | 0.855 → 0.433 | 0.1632 |
| C-grid, 23° | 4,993 | 7,153 | 2,967 | 19,232 | 0.0801 | 0.998 → 0.586 | 0.0112 |

The first C-grid trial allowed subdivisions based on orientation alone and produced a minimum metric mean ratio of 0.00197. The final implementation also requires a local metric quality check for splits. The final run rejected 235 splits on quality and maintained a minimum metric mean ratio of 0.0801. Low physical mean ratios remain, including near geometric constraints. Applying this mesh to a CFD case requires separate acceptance criteria for physical quality and near-wall resolution.

### Directional response and inflection

The median angular error of the fine metric axis in the peak-curvature band was 1.63°, 1.45° and 3.97°, respectively. The medians of $n^\top M_b n/(t^\top M_b t)$ in that band were 31.17, 23.93 and 10.24. Within 1.4 transition widths of the centreline, edge subsets within 25° of the normal or tangent gave the following physical length medians.

| Case | Near-normal edges, source → adapted | Near-tangential edges, source → adapted |
| --- | ---: | ---: |
| Square, 19° | 0.0674 → 0.0311 | 0.0625 → 0.0780 |
| Square, 67° | 0.0746 → 0.0269 | 0.0634 → 0.0767 |
| C-grid, 23° | 0.00401 → 0.201 | 0.00402 → 0.514 |

The squares show normal refinement and tangential coarsening for both rotations. The C-grid directional subsets include the densely sampled airfoil neighborhood. Both directional medians grew there, and adapted near-normal edges were about 2.6 times shorter than adapted near-tangential edges. The sample reflects redistribution of the initial geometric refinement as well as the prescribed transition.

The exact normal curvature at the tanh inflection is zero. The recovered, interpolated normal curvature at each specified centre was 0.137, 0.0103 and −0.131 in that case's density per length squared. Recovery averages across finite neighborhoods, and the nonuniform C-grid samples the transition asymmetrically. The reference pipeline retains the floor and resolved flank curvature without an added indicator.

### Geometry, tooling and remaining verification

The C-grid kept one connected fluid region with two boundary loops and Euler characteristic zero. Its area was 1827.22749075084 chord squared to the displayed precision. All 643 fixed source coordinates and 38 embedded wake edges were retained exactly. Marker segment counts changed from 524/62/20 to 524/120/40 through splits on straight portions of the farfield and outlet. All three marker lengths were unchanged to reported precision.

ParaView 5.13.2 loaded the triangle and line VTUs for all cases and returned the expected cell counts, fields and tensor dimensions. Gmsh 4.15.2 independently checked coherence of an MSH2 triangle snapshot made from the exported C-grid coordinates and connectivity. The SU2 ASCII round trips passed in the Python interface. Validation with a real SU2 restart export, loading the adapted mesh in the SU2 8.5.0 executable, and assessing its CFD performance remain work for Mission 03. The cloud timings do not establish native macOS performance.
