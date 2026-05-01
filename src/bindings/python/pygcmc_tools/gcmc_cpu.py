from __future__ import annotations

import os
import stat
import sys
from importlib import resources


def main() -> None:
    binary = resources.files("pygcmc_tools").joinpath("bin", "gcmc_cpu")
    binary_path = os.fspath(binary)
    current_mode = os.stat(binary_path).st_mode
    if not current_mode & stat.S_IXUSR:
        os.chmod(binary_path, current_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    os.execv(binary_path, [binary_path, *sys.argv[1:]])


if __name__ == "__main__":
    main()
