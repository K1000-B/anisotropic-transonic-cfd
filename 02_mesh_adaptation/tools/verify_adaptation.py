#!/usr/bin/env python3
"""Reproduce the manufactured and Mission 01 adaptation evidence."""

import argparse
import cProfile
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform
import pstats
import subprocess
from time import perf_counter

import numpy as np

from cfd_adapt import GeometryConstraints, MetricParameters, RemeshParameters, adapt
from cfd_adapt.diagnostics import distribution, summarise, write_diagnostics
from cfd_adapt.field import P1MetricField
from cfd_adapt.io import read_su2, write_su2
from cfd_adapt.mesh import edge_key
from cfd_adapt.metric import recover_hessian
from cfd_adapt.verification import grid, tanh_transition


ROOT = Path(__file__).resolve().parents[2]
BASELINE = ROOT / "01_baseline_mesh/reference/baseline.su2"


def baseline_constraints(mesh):
    """Protect the existing unit-chord wake embedding and curved CAD traces."""
    x, y = mesh.coordinates.T
    seam = np.flatnonzero((np.abs(y) < 1e-12) & (x >= 1) & (x <= 30))
    seam = seam[np.argsort(x[seam])]
    protected = frozenset(edge_key(a, b) for a, b in zip(seam[:-1], seam[1:]))
    if len(seam) < 3 or x[seam[0]] != 1 or x[seam[-1]] != 30 or not protected <= mesh.edges.keys():
        raise ValueError("expected the embedded wake of the reference C-grid")
    straight = set()
    for edge, marker in mesh.boundary.items():
        a, b = list(edge)
        if (marker == "outlet" or (marker == "farfield" and y[a] == y[b] and abs(y[a]) == 20)
                or (marker == "airfoil" and x[a] == x[b] == 1)):
            straight.add(edge)
    return GeometryConstraints(protected_edges=protected, splittable_boundary_edges=frozenset(straight))


def case_settings(name):
    if name == "baseline":
        mesh = read_su2(BASELINE)
        settings = MetricParameters(0.0005, 0.01, 4.0, 12, target_nodes=4200)
        return mesh, 23., np.array([1.7, 0.2]), 0.65, settings, baseline_constraints(mesh)
    mesh = grid(16, jitter=0.15)
    settings = MetricParameters(0.05, 0.023, 0.35, 7, target_nodes=280)
    angle = float(name.removeprefix("square-"))
    return mesh, angle, np.array([0.5, 0.5]), 0.09, settings, GeometryConstraints(splittable_boundary_edges=None)


def quadratic_verification():
    h = np.array([[1., 0.4], [0.4, -0.6]])
    evidence = []
    for jitter in (0., 0.17):
        for n in (12, 24, 48):
            mesh = grid(n, jitter)
            rho = 1 + 0.5 * np.einsum("ni,ij,nj->n", mesh.coordinates, h, mesh.coordinates)
            errors = np.linalg.norm(recover_hessian(mesh, rho).hessian - h, axis=(1, 2))
            boundary = np.unique(list(mesh.boundary))
            interior = np.all((mesh.coordinates > 2.1/n) & (mesh.coordinates < 1 - 2.1/n), axis=1)
            evidence.append({"cells_per_side": n, "jitter_fraction": jitter,
                             "area_weighted_rms": float(np.sqrt(np.average(errors ** 2, weights=mesh.nodal_areas))),
                             "boundary_mean": float(errors[boundary].mean()),
                             "interior_mean": float(errors[interior].mean())})
    return evidence


def direction_statistics(mesh, origin, normal, width):
    edges = np.array(list(mesh.edges))
    endpoints = mesh.coordinates[edges]
    vectors = endpoints[:, 1] - endpoints[:, 0]
    lengths = np.linalg.norm(vectors, axis=1)
    normal_fraction = abs(vectors @ normal) / lengths
    tangent_fraction = abs(vectors @ np.array([-normal[1], normal[0]])) / lengths
    in_band = abs((endpoints.mean(axis=1) - origin) @ normal) <= 1.4 * width
    return {
        label: {"count": len(subset), "length": distribution(subset) if len(subset) else None}
        for label, values in (("normal", normal_fraction), ("tangent", tangent_fraction))
        for subset in [lengths[in_band & (values >= np.cos(np.deg2rad(25)))]]
    }


def verify_geometry(mesh, result, constraints):
    output = result.mesh
    by_origin = {v: i for i, v in enumerate(result.source_vertex_indices) if v >= 0}
    rim = set(v for edge in mesh.boundary for v in edge)
    fixed = rim | set(constraints.fixed_vertices) | {v for edge in constraints.protected_edges for v in edge}
    assert all(v in by_origin and np.array_equal(mesh.coordinates[v], output.coordinates[by_origin[v]]) for v in fixed)
    assert all(edge_key(by_origin[a], by_origin[b]) in output.edges for a, b in constraints.protected_edges)
    # Boundary operations only split existing oriented line segments. Check
    # marker lengths independently and protect every non-splittable source edge.
    perimeter_errors = {}
    for name, edges in mesh.markers.items():
        def length(m, e):
            return np.linalg.norm(m.coordinates[e[:, 0]] - m.coordinates[e[:, 1]], axis=1).sum()
        perimeter_errors[name] = float(length(output, output.markers[name]) - length(mesh, edges))
        assert abs(perimeter_errors[name]) < 1e-12 * length(mesh, edges)
    if constraints.splittable_boundary_edges is not None:
        for (a, b), name in mesh.boundary.items():
            if (a, b) not in constraints.splittable_boundary_edges:
                assert output.boundary[edge_key(by_origin[a], by_origin[b])] == name
    return {"fixed_source_vertices": len(fixed), "protected_edges": len(constraints.protected_edges),
            "marker_length_changes": perimeter_errors}


def plots(path, source, result, target_nodes):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import EllipseCollection

    xy = source.coordinates
    m = result.metric
    d = result.diagnostics
    fig, axes = plt.subplots(3, 3, figsize=(13, 10), layout="constrained")

    def mesh_plot(axis, mesh, title):
        axis.triplot(*mesh.coordinates.T, mesh.triangles, lw=0.26, color="#233e4e")
        axis.set(title=title, aspect="equal", xlabel="x", ylabel="y")

    def field_plot(axis, values, title, label):
        plot = axis.tripcolor(*xy.T, source.triangles, values, shading="gouraud", cmap="viridis")
        fig.colorbar(plot, ax=axis, label=label, fraction=0.05)
        axis.set(title=title, aspect="equal", xlabel="x", ylabel="y")

    mesh_plot(axes[0, 0], source, "Source mesh")
    field_plot(axes[0, 1], result.density, "Manufactured density", "Density")
    field_plot(axes[0, 2], np.log10(m.curvature_eigenvalues[:, 1]), "Regularised curvature", "log₁₀ k₊")
    axis = axes[1, 0]
    axis.triplot(*xy.T, source.triangles, lw=0.2, color="#b2c1c9")
    # One node per spatial bin when possible keeps ellipses legible.
    span = np.ptp(xy, axis=0)
    bins = np.floor(18 * (xy - xy.min(axis=0)) / span).astype(int)
    _, indices = np.unique(bins, axis=0, return_index=True)
    orientation = np.rad2deg(np.arctan2(m.principal_directions[indices, 1, 0],
                                       m.principal_directions[indices, 0, 0]))
    ellipses = EllipseCollection(0.8 * m.target_sizes[indices, 0], 0.8 * m.target_sizes[indices, 1],
                                 orientation, units="xy", offsets=xy[indices], transOffset=axis.transData,
                                 facecolors="none", edgecolors="#b65134", linewidths=0.8)
    axis.add_collection(ellipses)
    axis.set(title="Metric ellipses (0.4× semiaxes)", aspect="equal", xlabel="x", ylabel="y")
    field_plot(axes[1, 1], np.log10(m.target_sizes[:, 1]), "Short principal size", "log₁₀ hₛ")
    mesh_plot(axes[1, 2], result.mesh, "Adapted mesh")
    a, b = d.initial_edge_lengths, d.final_edge_lengths
    bins = np.geomspace(min(a.min(), b.min()), max(a.max(), b.max()), 42)
    axes[2, 0].hist([a, b], bins=bins, density=True, label=["Source", "Adapted"], histtype="step", linewidth=1.4)
    axes[2, 0].axvline(1, color="#444", linestyle=":")
    axes[2, 0].set(xscale="log", xlabel="Metric edge length", ylabel="Probability density")
    axes[2, 0].legend()
    for metric_quality, before, after in (
        (False, d.initial_physical_quality, d.final_physical_quality),
        (True, d.initial_metric_quality, d.final_metric_quality),
    ):
        name = "Metric" if metric_quality else "Physical"
        axes[2, 1].hist([before, after], bins=np.linspace(0, 1, 24), density=True,
                        label=[f"{name}, source", f"{name}, adapted"], histtype="step", linewidth=1.2)
    axes[2, 1].set(xlabel="Triangle mean-ratio quality", ylabel="Probability density")
    axes[2, 1].legend(fontsize=8)
    counts = [len(source.coordinates), target_nodes, len(result.mesh.coordinates)]
    bars = axes[2, 2].bar(["Source", "Target", "Adapted"], counts, color=["#2e6d80", "#b1b9bc", "#b65134"])
    axes[2, 2].bar_label(bars, padding=3)
    axes[2, 2].set(ylabel="Vertices", ylim=(0, max(counts) * 1.15))
    fig.savefig(path / "verification.png", dpi=160)
    fig.savefig(path / "verification.pdf")
    plt.close(fig)


def run(name, output, level, profile, figures, passes):
    path = output / name
    path.mkdir()
    source, angle, origin, width, settings, constraints = case_settings(name)
    theta = np.deg2rad(angle)
    normal = np.array([np.cos(theta), np.sin(theta)])
    density, exact, distance = tanh_transition(source.coordinates, normal, origin, width)
    controls = RemeshParameters(max_passes=passes)
    start = perf_counter()
    if profile:
        profiler = cProfile.Profile()
        result = profiler.runcall(adapt, source, density, settings, controls, constraints)
        profiler.dump_stats(path / "profile.pstats")
        with (path / "profile.txt").open("w", encoding="utf-8") as stream:
            pstats.Stats(profiler, stream=stream).strip_dirs().sort_stats("cumulative").print_stats(30)
    else:
        result = adapt(source, density, settings, controls, constraints)
    wall = perf_counter() - start
    evidence = summarise(source, result, settings, controls)
    evidence["timings"]["adapt_wall_seconds"] = wall
    evidence["geometry"] = verify_geometry(source, result, constraints)
    error = np.linalg.norm(result.recovery.hessian - exact, axis=(1, 2))
    near_peak = (abs(distance) >= 0.40 * width) & (abs(distance) <= 0.90 * width)
    angle_error = np.rad2deg(np.arccos(np.clip(abs(result.metric.principal_directions[near_peak, :, 1] @ normal), 0, 1)))
    peak_metric = result.metric.bounded[near_peak]
    tangent = np.array([-normal[1], normal[0]])
    field = P1MetricField(source, result.metric.bounded)
    face, weights = field.locate(origin)
    recovered_at_inflection = weights @ np.einsum("i,nij,j->n", normal,
                                                  result.recovery.hessian[source.triangles[face]], normal)
    evidence["manufactured"] = {
        "angle_degrees": angle, "origin": origin.tolist(), "width": width, "base_density": 1.0, "jump": 0.6,
        "source_description": ("Mission 01 baseline.su2" if name == "baseline"
                               else "16 cells per side; 0.15 cell jitter; NumPy seed 41"),
        "source_sha256": hashlib.sha256(BASELINE.read_bytes()).hexdigest() if name == "baseline" else None,
        "area_weighted_hessian_rms": float(np.sqrt(np.average(error ** 2, weights=source.nodal_areas))),
        "peak_band_normal_direction_error_degrees": distribution(angle_error),
        "peak_band_normal_to_tangent_metric_eigenvalue_proxy": distribution(
            (normal @ peak_metric @ normal) / (tangent @ peak_metric @ tangent)
        ),
        "analytic_inflection_curvature": 0.0,
        "sampled_reconstructed_inflection_normal_curvature": float(recovered_at_inflection),
        "source_directional_edges": direction_statistics(source, origin, normal, width),
        "adapted_directional_edges": direction_statistics(result.mesh, origin, normal, width),
    }
    write_su2(path / "adapted.su2", result.mesh)
    reread = read_su2(path / "adapted.su2")
    assert np.array_equal(reread.coordinates, result.mesh.coordinates)
    assert np.array_equal(reread.triangles, result.mesh.triangles)
    assert reread.boundary == result.mesh.boundary
    evidence["adapted_su2_sha256"] = hashlib.sha256((path / "adapted.su2").read_bytes()).hexdigest()
    write_diagnostics(path, source, result, evidence, level=level)
    if figures:
        plots(path, source, result, settings.target_nodes)
    print(f"{name}: {len(source.coordinates)} -> {len(result.mesh.coordinates)} vertices, "
          f"{wall:.3f} s, {result.diagnostics.termination}", flush=True)
    return evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    cases = ["square-19", "square-67", "baseline"]
    parser.add_argument("output", type=Path)
    parser.add_argument("--case", choices=[*cases, "all"], default="all")
    parser.add_argument("--passes", type=int, default=8)
    parser.add_argument("--diagnostics", choices=["summary", "vtk", "full"], default="summary")
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--figures", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True)
    evidence = {
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "platform": platform.platform(), "profiled": args.profile,
                        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                        "worktree_dirty": bool(subprocess.check_output(
                            ["git", "status", "--porcelain"], cwd=ROOT, text=True).strip())},
        "quadratic_recovery": quadratic_verification(),
        "cases": {name: run(name, args.output, args.diagnostics, args.profile, args.figures, args.passes)
                  for name in (cases if args.case == "all" else [args.case])},
    }
    (args.output / "results.json").write_text(json.dumps(evidence, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
