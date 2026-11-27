# NACA 0012 baseline triangular mesh

Status: generated and geometrically checked, not CFD-validated. The user selected NACA 0012 on 2026-11-27. No flow conditions, shock position, turbulence model or wall-resolution requirement has been chosen.

## Geometry and provenance

For x/c = X in [0,1], the standard symmetric profile is defined by

`y/c = +/- 0.6*(0.2969*sqrt(X) - 0.1260*X - 0.3516*X^2 + 0.2843*X^3 - 0.1015*X^4)`.

Primary historical reference: Jacobs, Ward and Pinkerton, *The characteristics of 78 related airfoil sections from tests in the variable-density wind tunnel*, NACA-TR-460 (1933), [NTRS 19930091108](https://ntrs.nasa.gov/citations/19930091108). NTRS identifies it as a US Government work with public use permitted. The [NASA-origin Turbulence Modeling Resource](https://tmbwg.github.io/turbmodels/naca0012_val.html) explicitly gives the original formula and distinguishes its altered sharp-edge geometry.

We preserve the original finite trailing-edge thickness, `0.00252*c`, and close it with a straight vertical cap at X=1, split at its midpoint for the wake seam. We do not substitute the -0.1036 coefficient or the NASA extended/rescaled sharp-edge variant. This convention must be checked against the eventual experimental specimen before validation. No third-party coordinate file or protected document is redistributed: the analytical definition is evaluated directly in native Gmsh loops. No observed experimental data are generated here.

`geometry/naca0012_profile.geo` is reusable independently of the domain; its input/output variables are stated in its header. Cosine-spaced samples and interpolating splines approximate the analytical curves (301 points per side, shared leading edge). This is a documented geometric approximation, not an exact analytical CAD representation. The mesh boundary itself consists of straight segments.

## Domain and sizing

The C-domain has an upstream half ellipse centered at the leading edge, horizontal upper/lower extensions, and a vertical downstream boundary. Default upstream and vertical extents are equal, producing a semicircle. A shared embedded line joins the trailing-edge cap midpoint to the downstream midpoint. It is an internal conforming seam, never a physical boundary: this is an unstructured triangular C-layout, not a structured C-grid with index lines or separate duplicated cut faces.

All length controls below except Chord are relative to chord. Chord=1 is a geometric normalization, not a prescribed physical dimensional chord. AirfoilX=AirfoilY=0 locate the leading edge; incidence remains a later flow-setting decision.

| Parameter | Default | Purpose and rationale |
| --- | --- | --- |
| `Chord` | 1 | Unit chord normalization |
| `Upstream`, `Farfield` | 20, 20 | Initial substantial separation from the body; no domain-independence claim |
| `Downstream` | 30 | Outlet x-coordinate from the leading edge: 29 chords behind the trailing edge |
| `FarSize` | 2 | Coarse exterior; keeps cost low away from the geometric refinement |
| `AirfoilSize` | 0.004 | About 250 intervals per chord; also avoids poor triangles at the finite trailing-edge cap |
| `LeadingSize` | 0.002 | Extra nose resolution for strong curvature, independent of flow |
| `AirfoilInner`, `AirfoilOuter` | 0.03, 3 | Uniform near-body sizing followed by a broad linear transition |
| `LeadingInner`, `LeadingOuter` | 0.02, 0.15 | Local leading-edge refinement, returning to the other field outside this range |
| `GeometrySamples` | 301 | Spline support sampling, independent of mesh density |
| `DistanceSamples` | 1000 | Samples used by Gmsh's curve-distance field; not an exact distance oracle |
| `Algorithm` | 6 | Native Frontal-Delaunay triangulation; 5 is exposed as an alternative |

Two scalar `Distance`/`Threshold` pairs, combined by `Min`, control sizes. Point/curvature sizing and automatic boundary-size extension are disabled to avoid hidden controls. There is no anisotropic metric, boundary-layer extrusion, recombination, shock indicator or wake-size field. Ordinary triangles are not exactly equilateral, and the conforming seam constrains their edges; neither constitutes prescribed anisotropic refinement. Native Gmsh arrays, loops, conditions, expressions, splines, arcs/ellipses, embedding and physical groups implement all geometry and sizing. See the [official manual](https://gmsh.info/doc/texinfo/).

An initial 0.008 near-body size produced a 9.10-degree minimum angle around the finite-edge geometry. It was rejected during geometric development; the final uniform near-body size gives 22.31 degrees. This correction is not a mesh-convergence study. Domain dimensions and size settings remain initial engineering choices requiring later solution sensitivity checks.

## Reproduction

From the repository root, with the installed Gmsh on PATH:

```sh
mkdir -p runs/baseline/default
gmsh 01_baseline_mesh/geometry/naca0012_cgrid.geo -2 -format msh2 \
  -o runs/baseline/default/baseline.msh -nopopup
MPLCONFIGDIR=/private/tmp/cfd-mpl-baseline .PROJECT_CFDvenv/bin/python -B \
  01_baseline_mesh/tools/verify_mesh.py \
  runs/baseline/default/baseline.msh runs/baseline/default
gmsh runs/baseline/default/baseline.msh -check -nopopup
```

The verifier checks its interpreter before scientific imports, uses the Gmsh Python binding already bundled with the installed CLI, and writes SU2 through **Gmsh's native exporter**. It is not a mesh converter implementation. It also needs the existing NumPy/Matplotlib. The current binding lookup follows the installed Gmsh prefix (`bin/../lib/gmsh.py`); another installation layout would need explicit review. No dependencies were installed.

Open the main `.geo` in the Gmsh GUI to edit its parameter panel. For CLI overrides, append e.g. `-setnumber FarSize 1.5` or `-setnumber Upstream 15` to the generation command; no source rewrite is required. Non-default ellipse parsing was checked. The quantitative verifier and its plot limits currently intentionally qualify only the default unit-chord, zero-offset case; extend that limited analysis explicitly before using it for translated/scaled cases. The native `.geo` supports their generation.

## Results, metric definitions and acceptance

The reviewed `reference/` snapshot contains both mesh formats, `quality.json` and five figures. Total snapshot size is about 3.95 MB, with each mesh below 1 MB: a deliberate selected-reference exception to ordinary ignored `runs/` storage. Do not accumulate subsequent generated runs here. The JSON stores exact mesh SHA-256 values and Gmsh version. Logs and drafts remain in ignored run storage.

| Indicator | Observed |
| --- | --- |
| Exported nodes / linear triangles | 7154 / 13702 |
| `airfoil` / `farfield` / `outlet` segments | 524 / 62 / 20 |
| Edge length range / chord | 0.001223 to 2.597659 |
| Area range / chord squared | 1.018e-6 to 2.134431 |
| Mean-ratio quality minimum / median | 0.53837 / 0.99841 |
| Minimum angle | 22.31343 degrees |
| Maximum aspect ratio | 2.37346 |
| Maximum adjacent equivalent-size ratio | 1.97154 |

For each triangle of area A, side lengths a,b,c, perimeter P and inradius r=2A/P: quality is `4*sqrt(3)*A/(a*a+b*b+c*c)` (equilateral=1); aspect ratio is `max(a,b,c)/(2*sqrt(3)*r)` (equilateral=1). Angles use the cosine rule. Signed areas assume positive counterclockwise orientation. Size statistics use all triangle edge incidences, not unique edges; requested size fields are targets, not strict bounds on individual edge lengths. Adjacent-size jumps use `sqrt(max(A1,A2)/min(A1,A2))`. The JSON includes min, 1st percentile, median, 95th/99th percentiles and max.

Acceptance for this baseline uses positive areas, closed conforming boundaries, one connected component, no unused or duplicate exported nodes, only manifold edge incidences, quality above 0.5 and angles above 20 degrees. These are geometric screening criteria, not CFD accuracy guarantees. Duplicate screening rounds each coordinate to 12 decimals in the unit-chord case; it is not a general near-coincidence search. Euler characteristic V-E+T=0 corroborates one connected domain with one airfoil hole. Every one-sided edge is covered exactly once by a physical boundary marker; the seam has no extra one-sided boundary.

Boundary-node vertical error against the analytical profile is 1.59e-9 chord. The sampled analytical-to-upper-boundary-polyline distance is 3.16e-5 chord (20,001 cosine-distributed query points): this captures linear faceting, is a sampled estimate rather than a rigorous Hausdorff bound, and does not quantify the lower surface independently. No geometry-convergence study has been performed. Gmsh's native coherence check passed. Overall, profile, nose, trailing-edge and downstream images were visually inspected.

## SU2 interface and limitations

The export uses `NDIME=2`, only volume triangle type 5 and boundary line type 3. Export counts, coordinates and triangle connectivity were checked against MSH accounting for native node renumbering. Markers are exactly `airfoil`, `farfield`, `outlet`; the physical surface is `fluid`. Marker names do not prescribe solver boundary conditions; those remain to be justified in a later configuration. No internal-wake marker is exported.

The ASCII SU2 format is intended for the qualified SU2 8.5.0 environment. This is a structural/export compatibility check, **not** an actual solver load test: no SU2 configuration, driver initialization or CFD execution was performed. No convergence or physical validation claim follows from this mesh. The final experimental geometry, far-field independence and shock-resolution adequacy remain to be established.
