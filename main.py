#!/usr/bin/env python
"""
Wrapper for the champagne CLI.
Allows running the CLI out-of-the-box without first installing it with pip.
"""
import importlib.util
import os
import re
import sys
import tomllib

# add script directory to the path to allow champagne CLI to work out-of-the-box
# without the need to install it via pip first
SCRIPT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "src")
sys.path.append(SCRIPT_DIR)
from src.__main__ import main


def check_packages_installed(proj_file="pyproject.toml"):
    """Check whether required dependencies are installed.
    WARNING: Does not check for the required version.
    """

    with open(proj_file, "rb") as file:
        proj_data = tomllib.load(file)
    package_deps = {pkg.split()[0] for pkg in proj_data["project"]["dependencies"]}
    missing_pkgs = {pkg for pkg in package_deps if not importlib.util.find_spec(pkg)}
    if any(missing_pkgs):
        raise ImportError(
            f"Required Python packages are not installed: {', '.join(missing_pkgs)}"
        )


if (
    __name__ == "__main__"
):  # this block is adapted from the bin/champagne file created by pip
    check_packages_installed()
    sys.argv[0] = re.sub(r"(-script\.pyw|\.exe)?$", "", sys.argv[0])
    sys.exit(main())
