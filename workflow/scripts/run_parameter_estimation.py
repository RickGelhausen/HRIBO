#!/usr/bin/env python3
"""Serialize one DeepRibo cutoff publication and exec its R wrapper."""

import fcntl
import os
import shutil
import stat
import sys


R_STARTUP_OVERRIDES = (
    "R_ARCH",
    "R_DEFAULT_PACKAGES",
    "R_ENVIRON",
    "R_ENVIRON_USER",
    "R_HOME",
    "R_LIBS",
    "R_LIBS_SITE",
    "R_LIBS_USER",
    "R_PROFILE",
    "R_PROFILE_USER",
)


class LaunchError(RuntimeError):
    """Raised when the parameter-estimation process cannot be launched safely."""


def option_value(arguments, names, required=False):
    """Return the single value following any of *names*."""

    positions = [index for index, value in enumerate(arguments) if value in names]
    if len(positions) > 1:
        raise LaunchError("option supplied more than once: {}".format(names[-1]))
    if not positions:
        if required:
            raise LaunchError("required option is missing: {}".format(names[-1]))
        return None
    index = positions[0]
    if index + 1 >= len(arguments):
        raise LaunchError("option requires a value: {}".format(arguments[index]))
    return arguments[index + 1]


def output_directories(arguments):
    """Resolve every directory whose outputs belong to the published pair."""

    output = option_value(arguments, ("-o", "--out"), required=True)
    plot = option_value(arguments, ("-p", "--plot"))
    destination = option_value(arguments, ("-d", "--dest"))
    receipt = option_value(arguments, ("--receipt",))
    if plot is not None and destination is not None:
        raise LaunchError("use either --plot or --dest, not both")
    if plot is None:
        plot = (
            destination + ".png"
            if destination is not None
            else os.path.join(os.path.dirname(output), "s_curve.png")
        )

    directories = {
        os.path.realpath(
            os.path.abspath(os.path.expanduser(os.path.dirname(path) or "."))
        )
        for path in (output, plot, receipt)
        if path is not None
    }
    return sorted(directories)


def open_lock_directories(paths):
    """Open and exclusively lock output directories in a stable order."""

    descriptors = []
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    try:
        for path in paths:
            os.makedirs(path, exist_ok=True)
            descriptor = os.open(path, flags)
            mode = os.fstat(descriptor).st_mode
            if not stat.S_ISDIR(mode):
                raise LaunchError("output parent is not a directory: {}".format(path))
            descriptors.append(descriptor)
        for descriptor in descriptors:
            fcntl.flock(descriptor, fcntl.LOCK_EX)
            os.set_inheritable(descriptor, True)
    except Exception:
        for descriptor in descriptors:
            os.close(descriptor)
        raise
    return descriptors


def main(arguments):
    if not arguments:
        raise LaunchError(
            "usage: run_parameter_estimation.py WRAPPER.R [wrapper options]"
        )
    wrapper = os.path.realpath(os.path.expanduser(arguments[0]))
    wrapper_arguments = arguments[1:]
    if not os.path.isfile(wrapper):
        raise LaunchError("R wrapper does not exist: {}".format(arguments[0]))
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise LaunchError("Rscript is not available")

    descriptors = open_lock_directories(output_directories(wrapper_arguments))
    try:
        # The estimator and publication wrapper are checksum-/source-controlled.
        # --vanilla suppresses profile/Renviron files, but R still honors inherited
        # R_LIBS* and R_DEFAULT_PACKAGES values. Remove only those R startup and
        # package overrides while preserving the pinned container environment.
        environment = os.environ.copy()
        for name in R_STARTUP_OVERRIDES:
            environment.pop(name, None)
        os.execve(
            rscript,
            [rscript, "--vanilla", wrapper] + wrapper_arguments,
            environment,
        )
    finally:
        # Reached only when exec fails. Successful execution transfers these
        # inherited descriptors, and therefore the locks, to the R process.
        for descriptor in descriptors:
            os.close(descriptor)


if __name__ == "__main__":
    try:
        main(sys.argv[1:])
    except (LaunchError, OSError) as error:
        sys.stderr.write("run_parameter_estimation: error: {}\n".format(error))
        sys.exit(1)
