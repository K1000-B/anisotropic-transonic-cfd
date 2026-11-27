# Numerical verification

The [reproduction instructions](README.md) run the automated suite and the three manufactured-field cases. The [specification](mathematical-specification.md) maps individual checks to theory labels and implementation functions.

## Cases

- The square cases use a deterministic, slightly perturbed 16-by-16 cell grid and a tanh density transition through the centre. The two normal directions make angles of 19 and 67 degrees to the x-axis.
- The Mission 01 case uses the reference NACA 0012 C-grid with its original finite trailing edge and a tanh transition at 23 degrees through `(1.7, 0.2)` in chord units.
- The quadratic study compares reconstructed and exact Hessians on regular and jittered meshes, with separate interior and boundary statistics.

Verification checks each SU2 export by rereading coordinates, connectivity and markers. The C-grid check retains the polygonal airfoil, outer boundary and embedded wake. The manufactured inputs provide code-verification evidence. Flow-solution verification and physical validation belong to the later CFD stages.
