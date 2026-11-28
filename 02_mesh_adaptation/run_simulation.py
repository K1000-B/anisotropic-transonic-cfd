
"""
SU2 Docker Simulation Runner
============================

Execute a SU2 CFD simulation from the host operating system
using an existing Docker container.

Inputs:
    mesh_path   : SU2 mesh file
    config_path : SU2 configuration file

Execution procedure:

    1. Inspect the Docker container.
    2. Identify the shared directories.
    3. Convert host paths into container paths.
    4. Create a temporary SU2 configuration.
    5. Start the Docker container.
    6. Execute SU2_CFD inside the container.
    7. Wait for the simulation to finish.
    8. Stop the container if started by this program.
    9. Remove the temporary configuration.

The original mesh and configuration files are preserved.

No physical computation or post-processing is performed
outside SU2.

Only the Python standard library is required.
"""

import argparse
import json
import re
import subprocess
import tempfile

from pathlib import Path, PurePosixPath


# ============================================================
# 0. DOCKER CONFIGURATION
# ============================================================

# Name of the existing Docker container containing SU2.

CONTAINER_NAME = "su2-docker"

# Name of the SU2 executable inside the container.

SU2_EXECUTABLE = "SU2_CFD"


# ============================================================
# 1. INSPECT THE DOCKER CONTAINER
# ============================================================

def inspect_container(container_name):
    """
    Retrieve the Docker container configuration.

    Docker inspect provides information about:

        - Container state.
        - Shared volumes.
        - Host directory paths.
        - Container directory paths.
        - Container configuration.

    The output is returned as a Python dictionary.
    """

    result = subprocess.run(
        ["docker", "inspect", container_name],
        capture_output=True,
        text=True,
        check=True
    )

    # Docker inspect returns a JSON array containing
    # the configuration of the requested container.

    container_info = json.loads(result.stdout)[0]

    return container_info


# ============================================================
# 2. CONVERT HOST PATHS INTO CONTAINER PATHS
# ============================================================

def host_to_container_path(host_path, container_info):
    """
    Convert an absolute host path into its corresponding
    path inside the Docker container.

    Example
    -------
    Host directory:

        /Users/user/dockerData

    Container directory:

        /home/cfd/DATA

    Host file:

        /Users/user/dockerData/project/mesh.su2

    Container file:

        /home/cfd/DATA/project/mesh.su2

    The mounted directory is automatically identified using
    Docker inspect. No hard-coded volume path is required.
    """

    host_path = Path(host_path).resolve()

    matching_mounts = []

    # Retrieve the directories mounted inside Docker.

    for mount in container_info["Mounts"]:

        # Only host directory bind mounts can provide
        # direct access to files located on the Mac.

        if mount["Type"] != "bind":
            continue

        source = Path(mount["Source"]).resolve()

        destination = PurePosixPath(
            mount["Destination"]
        )

        # Determine whether the requested host file
        # belongs to this shared directory.

        try:

            relative_path = host_path.relative_to(source)

        except ValueError:

            # The requested path does not belong
            # to this particular mount.

            continue

        matching_mounts.append((
            source,
            destination,
            relative_path,
            mount.get("RW", False)
        ))

    # Docker cannot access arbitrary host files unless
    # they belong to a mounted directory.

    if not matching_mounts:

        raise ValueError(
            f"Path is not accessible inside Docker:\n"
            f"{host_path}"
        )

    # If several mounts contain the requested file,
    # select the most specific mounted directory.

    source, destination, relative_path, writable = max(
        matching_mounts,
        key=lambda mount: len(mount[0].parts)
    )

    if not writable:
        raise ValueError(f"Docker bind mount is read-only: {host_path}")

    # Reconstruct the path using the Linux path convention.

    container_path = destination.joinpath(
        *relative_path.parts
    )

    return str(container_path)


# ============================================================
# 3. PREPARE THE TEMPORARY SU2 CONFIGURATION
# ============================================================

def prepare_configuration(
    config_path,
    container_mesh_path
):
    """
    Read the original SU2 configuration and replace
    MESH_FILENAME with the path accessible inside Docker.

    All other configuration parameters remain unchanged.
    """

    configuration = Path(config_path).read_text(
        encoding="utf-8"
    )

    # Identify an active MESH_FILENAME definition.
    #
    # Commented configuration lines beginning with "%"
    # are not selected.

    pattern = (
        r"(?m)^[ \t]*MESH_FILENAME[ \t]*=[^\r\n]*"
    )

    configuration, replacements = re.subn(

        pattern,

        lambda match: (
            f"MESH_FILENAME= {container_mesh_path}"
        ),

        configuration
    )

    # Several active mesh definitions would make
    # the configuration ambiguous.

    if replacements > 1:

        raise ValueError(
            "Multiple MESH_FILENAME definitions found."
        )

    # Append the mesh definition if it is not already
    # present in the original configuration.

    if replacements == 0:

        configuration += (
            f"\nMESH_FILENAME= {container_mesh_path}\n"
        )

    return configuration


# ============================================================
# 4. START THE DOCKER CONTAINER
# ============================================================

def start_container(container_name):
    """
    Start an existing Docker container.

    This function does not create a new container and
    does not modify its configuration.
    """

    print(
        f"Starting Docker container: {container_name}",
        flush=True
    )

    subprocess.run(
        ["docker", "start", container_name],
        check=True
    )


# ============================================================
# 5. STOP THE DOCKER CONTAINER
# ============================================================

def stop_container(container_name):
    """
    Stop the Docker container without deleting it.

    The container, its installed software and its
    persistent data remain available for future simulations.
    """

    print(
        f"Stopping Docker container: {container_name}",
        flush=True
    )

    subprocess.run(
        ["docker", "stop", container_name],
        check=True
    )


# ============================================================
# 6. EXECUTE SU2 INSIDE DOCKER
# ============================================================

def execute_su2(
    container_name,
    container_config_path,
    container_working_directory,
    log_path=None
):
    """
    Execute SU2_CFD inside the running Docker container.

    The simulation is executed synchronously.

    Python waits until SU2_CFD finishes before continuing.

    Standard output and standard error are forwarded
    directly to the host terminal.
    """

    # Equivalent terminal command:
    #
    # docker exec \
    #     --workdir /home/cfd/DATA/project \
    #     su2-docker \
    #     SU2_CFD simulation.cfg

    command = [

        "docker",

        "exec",

        "-e",
        "PYTHONDONTWRITEBYTECODE=1",

        "--workdir",
        container_working_directory,

        container_name,

        "bash", "--noprofile", "--norc", "-c",
        # Activation is required by the documented non-interactive contract.
        # The unique PID marker identifies only this invocation on interruption.
        'source /home/cfd/.python_envs/python310.sh >/dev/null && '
        'printf "%s\\n" "$$" > "$1" && shift && exec "$@"',
        "bash",
        container_config_path + ".pid",

        SU2_EXECUTABLE,

        container_config_path
    ]

    print("\nStarting SU2 CFD simulation.\n")

    print(
        f"Container:     {container_name}"
    )

    print(
        f"Configuration: {container_config_path}"
    )

    print(
        f"Directory:     {container_working_directory}"
    )

    print(
        "\n" + "=" * 60,
        flush=True
    )

    # Stream to the terminal and optionally retain the same output for audit.
    # No shell is used on the host; paths remain separate argv entries.
    log_file = None
    process = None
    try:
        if log_path is not None:
            log_file = open(log_path, "x", encoding="utf-8")
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace", bufsize=1
        )
        for line in process.stdout:
            print(line, end="", flush=True)
            if log_file is not None:
                log_file.write(line)
                log_file.flush()
        returncode = process.wait()
        if returncode:
            raise subprocess.CalledProcessError(returncode, command)
    except BaseException:
        if process is not None:
            # docker exec client termination alone need not stop its solver.
            # Signal only the PID whose working directory and executable match
            # this invocation; never stop an already-running shared container.
            try:
                subprocess.run([
                    "docker", "exec", container_name, "bash", "-c",
                    'if read -r pid < "$1" && [[ "$pid" =~ ^[0-9]+$ ]] '
                    '&& [[ "$(readlink /proc/$pid/cwd)" == "$2" ]]; then '
                    'read -r -d "" program < /proc/$pid/cmdline; '
                    'if [[ "${program##*/}" == "$3" ]]; then '
                    'kill -TERM "$pid"; fi; fi',
                    "bash", container_config_path + ".pid",
                    container_working_directory, SU2_EXECUTABLE,
                ], check=False, timeout=15)
            except (OSError, subprocess.TimeoutExpired) as cleanup_error:
                print(f"Solver interruption cleanup uncertain: {cleanup_error}", flush=True)
            try:
                process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait()
                print("Check the retained run: solver shutdown is unconfirmed.", flush=True)
        raise
    finally:
        if process is not None and process.stdout is not None:
            process.stdout.close()
        if log_file is not None:
            log_file.close()

    print("\n" + "=" * 60)

    print(
        "\nSU2 execution completed successfully.",
        flush=True
    )


# ============================================================
# 7. COMPLETE SIMULATION WORKFLOW
# ============================================================

def run_simulation(mesh_path, config_path, *, log_path=None):
    """
    Run a complete SU2 simulation inside Docker.

    Parameters
    ----------
    mesh_path : str or Path
        Path to the SU2 mesh on the host operating system.

    config_path : str or Path
        Path to the SU2 configuration on the host
        operating system.

    log_path : str or Path, optional
        New host log file, streamed simultaneously to the terminal.

    Notes
    -----
    The original configuration remains unchanged.

    Simulation outputs are written into the directory
    containing the original configuration file, unless
    other output paths are specified in the configuration.

    The container is stopped automatically if it was
    started by this function.
    """

    # --------------------------------------------------------
    # 7.1. Validate host input paths
    # --------------------------------------------------------

    mesh_path = Path(mesh_path).resolve()

    config_path = Path(config_path).resolve()

    if not mesh_path.is_file():

        raise FileNotFoundError(
            f"Mesh not found: {mesh_path}"
        )

    if not config_path.is_file():

        raise FileNotFoundError(
            f"Configuration not found: {config_path}"
        )

    if mesh_path.suffix.lower() != ".su2":

        raise ValueError(
            "The input mesh must have a .su2 extension."
        )

    if config_path.suffix.lower() != ".cfg":

        raise ValueError(
            "The configuration must have a .cfg extension."
        )

    # --------------------------------------------------------
    # 7.2. Inspect the Docker container
    # --------------------------------------------------------

    container_info = inspect_container(
        CONTAINER_NAME
    )

    # Record whether the container was already running.
    #
    # A container that was already running may contain
    # another process that must not be interrupted.

    was_running = container_info["State"]["Running"]

    # --------------------------------------------------------
    # 7.3. Convert the host paths
    # --------------------------------------------------------

    container_mesh_path = host_to_container_path(

        mesh_path,

        container_info
    )

    # Execute the simulation from the original
    # configuration directory.
    #
    # Relative paths in the SU2 configuration therefore
    # refer to this directory.

    container_working_directory = host_to_container_path(

        config_path.parent,

        container_info
    )

    # --------------------------------------------------------
    # 7.4. Prepare the updated configuration
    # --------------------------------------------------------

    configuration = prepare_configuration(

        config_path,

        # Workflow snapshots use a safe basename even when their parent path
        # contains spaces. SU2 config values are not shell-quoted arguments.
        mesh_path.name if mesh_path.parent == config_path.parent
        and not any(character.isspace() for character in mesh_path.name)
        else container_mesh_path
    )

    temporary_config = None

    started_here = False

    try:

        # ----------------------------------------------------
        # 7.5. Create a temporary configuration on the host
        # ----------------------------------------------------

        # The temporary configuration is created next to
        # the original configuration.
        #
        # Because this directory is shared with Docker,
        # the temporary file becomes immediately accessible
        # inside the container.

        with tempfile.NamedTemporaryFile(

            mode="w",

            suffix=".cfg",

            prefix="simulation_",

            dir=config_path.parent,

            encoding="utf-8",

            delete=False

        ) as file:

            temporary_config = Path(file.name)

            file.write(configuration)

        # Convert the temporary configuration path
        # into its equivalent Docker path.

        container_config_path = host_to_container_path(

            temporary_config,

            container_info
        )

        # ----------------------------------------------------
        # 7.6. Start the container
        # ----------------------------------------------------

        if not was_running:

            start_container(
                CONTAINER_NAME
            )

            started_here = True

        else:

            print(
                f"Container {CONTAINER_NAME} "
                f"is already running.",
                flush=True
            )

        # ----------------------------------------------------
        # 7.7. Execute the CFD simulation
        # ----------------------------------------------------

        execute_su2(

            CONTAINER_NAME,

            container_config_path,

            container_working_directory,
            log_path=log_path
        )

    finally:

        # ----------------------------------------------------
        # 7.8. Stop Docker and clean temporary files
        # ----------------------------------------------------

        # The finally block is executed even if SU2 fails
        # or the simulation is interrupted.

        try:

            if started_here:

                stop_container(
                    CONTAINER_NAME
                )

        finally:

            # Always remove the temporary configuration.
            #
            # The original configuration and CFD output
            # files remain unchanged by this cleanup.

            if temporary_config is not None:

                temporary_config.unlink(
                    missing_ok=True
                )
                Path(str(temporary_config) + ".pid").unlink(missing_ok=True)

    # --------------------------------------------------------
    # 7.9. Print the final execution summary
    # --------------------------------------------------------

    print("\nSimulation workflow completed.")

    print(
        f"Output directory: {config_path.parent}"
    )


# ============================================================
# 8. COMMAND-LINE INTERFACE
# ============================================================

def main():
    """
    Execute the simulation directly from macOS.

    Usage:

        python run_simulation.py mesh.su2 config.cfg
    """

    parser = argparse.ArgumentParser(

        description=(
            "Run a SU2 CFD simulation inside Docker."
        )
    )

    parser.add_argument(

        "mesh_path",

        help="Path to the SU2 mesh file."
    )

    parser.add_argument(

        "config_path",

        help="Path to the SU2 configuration file."
    )

    args = parser.parse_args()

    try:

        run_simulation(

            args.mesh_path,

            args.config_path
        )

    except (
        OSError,
        ValueError,
        subprocess.CalledProcessError
    ) as error:

        parser.exit(
            1,
            f"Error: {error}\n"
        )


if __name__ == "__main__":

    main()
