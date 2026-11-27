// Standard NACA 0012 with its original finite trailing-edge thickness.
// Coordinates are evaluated natively; see README for provenance and conventions.
SetFactory("Built-in");
DefineConstant[
  Chord = {1, Min 1e-6, Name "Geometry/Chord"},
  AirfoilX = {0, Name "Geometry/Leading edge x"},
  AirfoilY = {0, Name "Geometry/Leading edge y"},
  Upstream = {20, Min 2, Name "Domain/Upstream chords"},
  Downstream = {30, Min 2, Name "Domain/Downstream x from leading edge in chords"},
  Farfield = {20, Min 2, Name "Domain/Vertical extent in chords"},
  FarSize = {2, Min 1e-6, Name "Mesh/Far field size over chord"},
  AirfoilSize = {0.004, Min 1e-6, Name "Mesh/Airfoil size over chord"},
  LeadingSize = {0.002, Min 1e-6, Name "Mesh/Leading edge size over chord"},
  AirfoilInner = {0.03, Min 0, Name "Mesh/Airfoil constant size distance"},
  AirfoilOuter = {3, Min 1e-6, Name "Mesh/Airfoil transition end distance"},
  LeadingInner = {0.02, Min 0, Name "Mesh/Leading edge constant size distance"},
  LeadingOuter = {0.15, Min 1e-6, Name "Mesh/Leading edge transition end distance"},
  GeometrySamples = {301, Min 21, Max 2001, Step 2, Name "Geometry/Samples per side"},
  DistanceSamples = {1000, Min 100, Step 100, Name "Mesh/Curve distance samples"},
  Algorithm = {6, Choices{5,6}, Name "Mesh/Triangulation algorithm"}
];
If(Chord <= 0 || LeadingSize <= 0 || LeadingSize > AirfoilSize || AirfoilSize > FarSize || AirfoilInner < 0 || LeadingInner < 0 || AirfoilInner >= AirfoilOuter || LeadingInner >= LeadingOuter || Upstream <= 0 || Farfield <= 0 || Downstream <= 1 || GeometrySamples < 21 || GeometrySamples != Floor(GeometrySamples) || DistanceSamples < 2 || DistanceSamples != Floor(DistanceSamples))
  Error("Invalid scale or size transition parameters"); Abort;
EndIf
Include "naca0012_profile.geo";
center = newp; Point(center) = {AirfoilX,AirfoilY,0};
top = newp; Point(top) = {AirfoilX,AirfoilY+Farfield*Chord,0};
left = newp; Point(left) = {AirfoilX-Upstream*Chord,AirfoilY,0};
bottom = newp; Point(bottom) = {AirfoilX,AirfoilY-Farfield*Chord,0};
rightU = newp; Point(rightU) = {AirfoilX+Downstream*Chord,AirfoilY+Farfield*Chord,0};
rightM = newp; Point(rightM) = {AirfoilX+Downstream*Chord,AirfoilY,0};
rightL = newp; Point(rightL) = {AirfoilX+Downstream*Chord,AirfoilY-Farfield*Chord,0};
arcU = newl; arcL = arcU+1;
If(Abs(Upstream-Farfield) < 1e-12)
  Circle(arcU) = {top,center,left}; Circle(arcL) = {left,center,bottom};
Else
  If(Upstream > Farfield)
    major = left;
  Else
    major = top;
  EndIf
  Ellipse(arcU) = {top,center,major,left}; Ellipse(arcL) = {left,center,major,bottom};
EndIf
wallU = newl; Line(wallU) = {rightU,top};
wallL = newl; Line(wallL) = {bottom,rightL};
outL = newl; Line(outL) = {rightL,rightM};
outU = newl; Line(outU) = {rightM,rightU};
outer = newll; Curve Loop(outer) = {wallU,arcU,arcL,wallL,outL,outU};
inner = newll; Curve Loop(inner) = {u,capU,capL,-l};
fluid = news; Plane Surface(fluid) = {outer,inner};
// A single shared embedded seam: no duplicate wake nodes or physical wake marker.
wake = newl; Line(wake) = {te,rightM}; Curve{wake} In Surface{fluid};
Field[1] = Distance; Field[1].CurvesList = {u,l,capU,capL}; Field[1].Sampling = DistanceSamples;
Field[2] = Threshold; Field[2].InField = 1;
Field[2].SizeMin = AirfoilSize*Chord; Field[2].SizeMax = FarSize*Chord;
Field[2].DistMin = AirfoilInner*Chord; Field[2].DistMax = AirfoilOuter*Chord;
Field[3] = Distance; Field[3].PointsList = {nose};
Field[4] = Threshold; Field[4].InField = 3;
Field[4].SizeMin = LeadingSize*Chord; Field[4].SizeMax = FarSize*Chord;
Field[4].DistMin = LeadingInner*Chord; Field[4].DistMax = LeadingOuter*Chord;
Field[5] = Min; Field[5].FieldsList = {2,4}; Background Field = 5;
Mesh.MeshSizeFromPoints = 0; Mesh.MeshSizeFromCurvature = 0; Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeMin = LeadingSize*Chord; Mesh.MeshSizeMax = FarSize*Chord;
Mesh.Algorithm = Algorithm; Mesh.ElementOrder = 1; Mesh.RecombineAll = 0;
Mesh.MshFileVersion = 2.2; Mesh.Binary = 0; Mesh.SaveAll = 0;
Physical Curve("airfoil") = {u,capU,capL,l};
Physical Curve("farfield") = {wallU,arcU,arcL,wallL};
Physical Curve("outlet") = {outL,outU};
Physical Surface("fluid") = {fluid};
