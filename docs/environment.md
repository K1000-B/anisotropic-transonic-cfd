# Software environment qualification

Historical qualification date: 2026-11-27. Those environment diagnostics performed no dependency installation, image rebuild, container configuration change, mesh generation or CFD simulation. Baseline mesh generation was subsequently qualified in its own stage. Mission 03 (2026-11-28) uses documentation/source inspection only: no interpreter, dependency, Docker or SU2 execution was checked anew.

## Intended workflow environment (implementation only)

The new host entry point is `03_adaptive_workflow/adaptive_workflow.py`. It and the elementary extraction/adaptation programs execute natively on macOS using this checkout's existing `.PROJECT_CFDvenv`. NumPy is required for extraction/adaptation/orchestration and Matplotlib for non-interactive mesh figures; the runner otherwise uses the Python standard library plus the existing Docker CLI. These are actual direct requirements, not a claim that the complete package inventory is a portable requirements specification. No additional package was introduced or installed. Dependency availability and runtime compatibility were not reverified during Mission 03.

Only CFD resolution executes inside the existing `su2-docker` container, with expected SU2 8.5.0. Cases/results require an existing writable host bind mount accessible to the container. The runner discovers mounts at future runtime rather than embedding host paths, selects the most specific match and rejects read-only access. External original inputs/dependencies are copied into case/run/solve snapshots under the mounted project; host numerical processing does not use container Python, and no host dependency repairs the high-level container `SU2` package.

The workflow adds fresh-start local configurations, ASCII restart/VTU/history checks, plots, PVD collections and timing/history/provenance. Implementation does not establish Docker launch compatibility, installed SU2 output semantics, format interoperability or numerical correctness. See the [workflow guide](../03_adaptive_workflow/README.md) for the separate future execution/testing procedure. Its commands must not be run as part of this implementation mission.

## Versions and execution domains

| Component | Observed version | Domain |
| --- | --- | --- |
| Host project Python | CPython 3.14.7 | macOS ARM64, `.PROJECT_CFDvenv/` |
| Gmsh | `4.15.2-git` (Homebrew package 4.15.2) | macOS ARM64 CLI |
| ParaView / pvpython / pvbatch | 6.1.1 | macOS ARM64 application bundle |
| ParaView bundled Python | CPython 3.12.7 | Application interpreter |
| Docker CLI / Engine | 29.8.0 / 29.8.0 | macOS ARM64 client / Linux ARM64 engine |
| Docker Desktop | 4.91.0 (239619) | Local runtime |
| SU2 | 8.5.0 `Harrier` | Linux AMD64 container, under emulation |
| SU2 wrapper Python | CPython 3.10.19 | Container venv, x86_64 |
| Container default `python3` | CPython 3.13.3 | System Python, not wrapper Python |
| Open MPI / mpi4py | 4.1.5 / 4.1.2 | Container only |
| Git / GitHub CLI | 2.55.0 / 2.102.0 | Host; authenticated account `K1000-B` |

## Host Python

The user-designated environment is `.PROJECT_CFDvenv/`, with system site-packages disabled. Interpreter identity and prefix were checked before execution. `pip check` reported no broken requirements. Imports of NumPy 2.5.5, Matplotlib 3.11.2 and Trimesh 5.2.2 passed; no scientific computation or figure was produced.

[Host package inventory](host-python-inventory.txt) records all 13 installed distributions with exact versions. It is an observed baseline, not a declaration that all packages are required, a portable lock file or complete wheel/build provenance. Define actual direct dependencies and a suitable lock when real code requires them. Check Python 3.14 support for each future dependency; compatibility of absent libraries is unproven.

From the repository root, verify identity before project execution or installation:

```sh
.PROJECT_CFDvenv/bin/python -B - <<'PY'
import pathlib
import sys
expected = pathlib.Path('.PROJECT_CFDvenv').absolute()
assert pathlib.Path(sys.prefix) == expected
assert pathlib.Path(sys.executable).parent == expected / 'bin'
print(sys.executable, sys.version)
PY
```

Activation with `source .PROJECT_CFDvenv/bin/activate` is optional; explicit interpreter paths are preferred in automation. Use only this interpreter's `-m pip` for justified host dependency changes. No additional dependency is needed now.

## Docker identity and resources

- Container: `su2-docker`; user/group `cfd`, UID/GID `1001:1001`; working directory `/home/cfd`; command `/bin/bash`, no configured entrypoint. Initially stopped, started for authorized diagnostics and left running.
- Image tag: `dockerappsdaer/su2_pywrapper:v1`.
- Image ID: `sha256:9b029bab47b6f929cc3a7778bdc996449fccb111494ee4f965489b32527c7616`.
- Repository digest: `dockerappsdaer/su2_pywrapper@sha256:9b029bab47b6f929cc3a7778bdc996449fccb111494ee4f965489b32527c7616`.
- Image: `linux/amd64`; engine: Linux `aarch64`; host: macOS `arm64`. Container execution reports `x86_64` and works under emulation. Native timing and scaling are unverified.
- Host: 10 logical/physical CPUs reported, 16 GiB RAM. Docker VM: 4 CPUs and 5,155,713,024 bytes RAM (about 4.80 GiB). Container `nproc` reports 4; no extra CPU quota, CPU set or memory cap. Shared memory: 64 MiB. The VM budget is shared with other containers and is not an approved simulation allocation.

The container has pre-existing writable-layer files. An image digest alone does not describe runtime state. No pre-existing case or result was imported, modified or interpreted as a project result. Observed SHA-256 values:

| Artifact | SHA-256 |
| --- | --- |
| `SU2_CFD` | `14a5f6cbaf5f480447e4370e96cdccfbedaaf5bc1806e12ca5ea881238dc2c77` |
| `_pysu2.so` | `ecfa24bd0bdb548d1ac481d0ced964ed8bccb124448c8fce5e6cbe944945c420` |

Exact source revision, build flags and OpenMP capability were not established. The banner does not prove an unmodified upstream build.

## Mounts and permissions

The existing read/write bind mount maps the host parent directory of this checkout to `/home/cfd/DATA`; this project is visible as `/home/cfd/DATA/PROJECT_CFD`. No mount change was needed. Sibling host directories are exposed too; automation must confine writes to the selected project run directory.

Machine-specific absolute paths and application locations are recorded only in ignored `local/machine-notes.md` and `local/tool-paths.env`. Public commands use relative host paths or local variables. Container paths here describe this image, not private host locations.

The project and README appear as `1001:1001` inside Docker, with directory mode 755 and README mode 644. Host/container README checksums matched. A host-created temporary marker was read in Docker and a container-created marker was read on the host; both and their temporary directory were removed. This proves current project read/write access, not every future file's permissions. No ownership changes were made.

Do not execute the macOS venv inside Linux or load Linux wrapper libraries into host Python. Exchange files between runtimes. Actual mesh, solution and restart compatibility remains untested: no format round trip was performed.

## SU2, MPI and wrapper

`SU2_RUN=/home/cfd/Code/SU2_PYWRAPPER/bin` and `SU2_HOME=/home/cfd/Code/SU2_PYWRAPPER` are configured. Binaries `SU2_CFD`, `SU2_DEF`, `SU2_DOT`, `SU2_GEO` and `SU2_SOL` exist in `SU2_RUN`; auxiliary Python scripts are present but not all qualified. `SU2_CFD --help` returned version 8.5.0 without reading a case or solving. Linked MPI libraries resolved. Do not assume optional features or additional executable variants from the version alone.

Open MPI binaries are in `/home/cfd/Code/OMPI/bin`, libraries in `/home/cfd/Code/OMPI/lib`. Two-rank `hostname` and mpi4py collective `allgather` diagnostics passed. No MPI CFD or scaling test was performed. ParaView's host MPI and the container MPI are separate installations.

Raw `docker exec` can locate SU2 and display its help, but does not run interactive `.bashrc` activation. System `python3` is 3.13.3; plain `python` was absent from that raw PATH. Source `/home/cfd/.python_envs/python310.sh` explicitly for the wrapper: it selects `/home/cfd/venvs/su2-py310/bin/python`, sets `PYTHONPATH` to `SU2_RUN`, and supplies the Python 3.10 and MPI library paths.

After activation, `import pysu2` succeeded with its compiled extension and `CSinglezoneDriver` was available; no driver was instantiated. Wrapper-venv distributions: mpi4py 4.1.2, packaging 26.3, pip 26.2.1, setuptools 84.0.0, wheel 0.48.0. **High-level `import SU2` failed: `ModuleNotFoundError: numpy`.** High-level utilities remain unqualified until a justified, explicitly approved container change. No dependency was added. Low-level wrapper-driven solving is also untested.

## Non-interactive launch contract

This activation pattern was historically smoke-tested with diagnostics, before the current orchestration existed. The new runner implements the same explicit activation with additional output streaming and invocation PID handling; subsequent user-run evidence is recorded in [project status](project-status.md). The following historical diagnostic command is retained as reference, not a Mission 03 instruction to execute it. From the repository root:

```sh
. local/tool-paths.env
# Inspect state first; if stopped, docker start su2-docker is authorized.
docker inspect --format '{{.State.Running}}' su2-docker
docker exec -e PYTHONDONTWRITEBYTECODE=1 \
  --workdir "$CFD_CONTAINER_PROJECT" su2-docker \
  bash --noprofile --norc -c \
  'source /home/cfd/.python_envs/python310.sh >/dev/null && exec "$@"' \
  bash SU2_CFD --help
```

This exact pattern passed both SU2 help and `python -B -c 'import pysu2'`. It needs no TTY, forwards arguments separately and returns the executed command's status. [Docker exec documentation](https://docs.docker.com/reference/cli/docker/container/exec/) describes its running-container requirement, inherited environment and working-directory option. The ignored local configuration must be recreated on another machine from its inspected mounts; it is intentionally not a public dependency.

The implemented runner discovers mounted paths through container inspection, selects the solve directory with `--workdir`, sources the documented activation script and passes `SU2_CFD` and its temporary configuration as separate arguments. It runs serially; no MPI CFD launcher is added. Host orchestration retains configs/manifests, and optional runner logs stream simultaneously to the terminal. A unique PID marker supports a guarded interruption-cleanup attempt without stopping an initially running container. Full integration and interruption behavior remain unqualified; one user-run baseline does not qualify every lifecycle branch. Do not use interactive aliases or install host packages to repair container imports.

Verify every future option against the [SU2 8.5.0 template](https://github.com/su2code/SU2/blob/v8.5.0/config_template.cfg) and actual build when needed. The observed help agrees with the [8.5.0 main program](https://github.com/su2code/SU2/blob/v8.5.0/SU2_CFD/src/SU2_CFD.cpp); no numerical setting is prescribed.

## Gmsh and ParaView automation

Gmsh version, build-info and help checks passed. It reports `4.15.2-git`, MacOSARM, Parser, OpenCASCADE 7.9.3, OpenMP, Bamg and quadrilateral-meshing features among its build options. These features do not change the project's own adaptation scope.

Future workflows will invoke native parametric `.geo` programs with the verified CLI options `-2`, `-format`, `-o`, `-setnumber` and `-nopopup` as needed. `su2` is listed as an export format. Real geometry will define the parameter names. No geometry parse, mesh generation or export was attempted. Consult the [official Gmsh manual](https://gmsh.info/doc/texinfo/).

`paraview` is on PATH. `pvpython` and `pvbatch` are absent from PATH but present in the application's `Contents/bin`; both returned 6.1.1 and ran a temporary script importing `paraview.simple`, reporting Python 3.12.7 on ARM64. No pipeline, rendering or dataset loading was requested; `--dr` skipped user application settings. Use these executables via local `PARAVIEW_BIN`, not project-venv imports. See the [official batch scripting guide](https://docs.paraview.org/en/v6.1.1/Tutorials/SelfDirectedTutorial/batchPythonScripting.html).

## Remaining limitations

Runtime basics and mount interoperability passed; production CFD, restart handling, format interoperability, numerical correctness, ParaView rendering and parallel visualization remain untested. Python 3.14 support of future dependencies is unproven. Emulation, VM resources and shared memory may constrain future cases; no resource change or rebuild is warranted yet. Determine whether high-level SU2 tools are actually needed before proposing a container change.

Docker socket and ParaView MPI/socket checks required execution outside the sandbox. Automatic approval initially rejected Docker start; explicit user confirmation resolved it. Future sessions must respect their own permissions. Continuity remains repository files and Git; no global Codex setting, custom skill, hook or scheduled task was changed.
