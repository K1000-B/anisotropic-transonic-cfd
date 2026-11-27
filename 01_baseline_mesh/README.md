# Baseline mesh

Own the baseline geometry definition, native parametric Gmsh `.geo` programs, boundary/physical-group definitions, meshing parameters and mesh-quality documentation. No geometry or mesh exists yet.

Future inputs are approved airfoil data, domain choices and explicit mesh parameters. Future outputs are baseline triangular meshes, boundary identifiers and quality summaries consumed by the adaptation and CFD workflows. Generate large meshes under root `runs/`, preserving their inputs and provenance; review any reference mesh before tracking it.

Keep geometry reusable and separate from adjustable mesh controls. Later cases must reference the same definitions rather than copy them. Use Gmsh's language for its native operations; do not enumerate geometry through Python unnecessarily.
