"""Command-line adapter for a single nodal density field."""

import argparse
import json
from pathlib import Path

from . import MetricParameters, RemeshParameters, adapt
from .diagnostics import summarise, write_diagnostics
from .io import constraints_from_ids, read_density_csv, read_su2, write_su2


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path, help="source SU2 triangular mesh")
    parser.add_argument("solution", type=Path, help="nodal density CSV with PointID and/or x,y")
    parser.add_argument("config", type=Path, help="JSON configuration")
    parser.add_argument("output", type=Path, help="new output directory")
    args = parser.parse_args(argv)
    try:
        config = json.loads(args.config.read_text(encoding="utf-8"))
        allowed = {"metric", "remeshing", "constraints", "diagnostics", "solution"}
        if unexpected := config.keys() - allowed:
            raise ValueError(f"unknown configuration keys {sorted(unexpected)}")
        source = read_su2(args.mesh)
        density = read_density_csv(args.solution, source, **config.get("solution", {}))
        metric = MetricParameters(**config["metric"])
        remeshing = RemeshParameters(**config.get("remeshing", {}))
        constraints = constraints_from_ids(source, **config.get("constraints", {}))
        level = config.get("diagnostics", "summary")
        if level not in {"summary", "vtk", "full"}:
            raise ValueError("diagnostic level must be summary, vtk or full")
        if args.output.exists():
            raise ValueError("output directory already exists")
        result = adapt(source, density, metric, remeshing, constraints)
        summary = summarise(source, result, metric, remeshing)
        args.output.mkdir(parents=True)
        write_su2(args.output / "adapted.su2", result.mesh)
        write_diagnostics(args.output, source, result, summary, level=level)
        (args.output / "config.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    except (OSError, ValueError, TypeError, KeyError, ArithmeticError) as error:
        parser.exit(2, f"cfd-adapt: {error}\n")
    print(f"{len(source.coordinates)} -> {len(result.mesh.coordinates)} vertices, "
          f"complexity {result.metric.bounded_complexity:.10g}, {result.diagnostics.termination}")


if __name__ == "__main__":
    main()
