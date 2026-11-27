// Reusable analytical profile: requires Chord, AirfoilX, AirfoilY, GeometrySamples.
// Outputs nose, te, u, l, capU, capL, upper[], lower[].
// Native cosine spacing resolves the nose geometry independently of mesh size.
upper[] = {}; lower[] = {};
For i In {0:GeometrySamples-1}
  xi = (1-Cos(Pi*i/(GeometrySamples-1)))/2;
  yt = 0.6*(0.2969*Sqrt(xi)-0.1260*xi-0.3516*xi^2+0.2843*xi^3-0.1015*xi^4);
  p = newp; Point(p) = {AirfoilX+Chord*xi, AirfoilY+Chord*yt, 0};
  upper[] += {p};
  If(i == 0)
    lower[] += {p};
  Else
    p = newp; Point(p) = {AirfoilX+Chord*xi, AirfoilY-Chord*yt, 0};
    lower[] += {p};
  EndIf
EndFor
nose = upper[0]; teU = upper[GeometrySamples-1]; teL = lower[GeometrySamples-1];
u = newl; Spline(u) = {upper[]};
l = newl; Spline(l) = {lower[]};
te = newp; Point(te) = {AirfoilX+Chord, AirfoilY, 0};
capU = newl; Line(capU) = {teU,te};
capL = newl; Line(capL) = {te,teL};
