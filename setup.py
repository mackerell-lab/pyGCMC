from __future__ import annotations

import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext


PROJECT_ROOT = Path(__file__).resolve().parent


class CMakeExtension(Extension):
    def __init__(self, name: str, source_dir: str = ".") -> None:
        super().__init__(name, sources=[])
        self.source_dir = str((PROJECT_ROOT / source_dir).resolve())


class CMakeBuild(build_ext):
    def build_extension(self, ext: CMakeExtension) -> None:
        build_temp = Path(self.build_temp).resolve()
        build_temp.mkdir(parents=True, exist_ok=True)

        extension_dir = Path(self.get_ext_fullpath(ext.name)).resolve().parent
        extension_dir.mkdir(parents=True, exist_ok=True)

        cmake_args = [
            f"-DCMAKE_BUILD_TYPE={self._build_type()}",
            f"-DPython_EXECUTABLE={sys.executable}",
            f"-DPython3_EXECUTABLE={sys.executable}",
        ]
        cmake_args.extend(self._pybind11_args())
        cmake_args.extend(self._split_env_args("CMAKE_ARGS"))

        build_args = ["--config", self._build_type()]
        build_args.extend(self._split_env_args("CMAKE_BUILD_ARGS"))
        if not any(arg.startswith("-j") or arg.startswith("--parallel") for arg in build_args):
            build_args.extend(["--parallel", os.environ.get("CMAKE_BUILD_PARALLEL_LEVEL", "4")])

        subprocess.check_call(["cmake", ext.source_dir, *cmake_args], cwd=build_temp)
        subprocess.check_call(
            ["cmake", "--build", ".", "--target", "pygcmc", "gcmc_cpu", *build_args],
            cwd=build_temp,
        )

        self._copy_extension(build_temp, extension_dir, ext)
        self._copy_cli(build_temp)

    def _build_type(self) -> str:
        if self.debug or os.environ.get("DEBUG") == "1":
            return "Debug"
        return os.environ.get("CMAKE_BUILD_TYPE", "Release")

    def _pybind11_args(self) -> list[str]:
        try:
            import pybind11
        except ImportError:
            return []
        cmake_dir = Path(pybind11.get_cmake_dir())
        if cmake_dir.exists():
            return [f"-Dpybind11_DIR={cmake_dir}"]
        return []

    def _split_env_args(self, name: str) -> list[str]:
        value = os.environ.get(name, "").strip()
        return value.split() if value else []

    def _copy_extension(self, build_temp: Path, extension_dir: Path, ext: CMakeExtension) -> None:
        suffix = self.get_ext_filename(ext.name).removeprefix(ext.name)
        module_dir = build_temp / "modules" / "bindings"
        expected = module_dir / f"{ext.name}{suffix}"
        if expected.exists():
            source = expected
        else:
            candidates = sorted(module_dir.glob(f"{ext.name}*.so"))
            if not candidates:
                raise RuntimeError(f"CMake did not produce the {ext.name} extension")
            source = candidates[0]
        if not source.exists():
            raise RuntimeError(f"CMake did not produce the {ext.name} extension")
        shutil.copy2(source, extension_dir / f"{ext.name}{suffix}")

    def _copy_cli(self, build_temp: Path) -> None:
        source = build_temp / "bin" / "gcmc_cpu"
        if not source.exists():
            raise RuntimeError("CMake did not produce the gcmc_cpu executable")
        package_bin = Path(self.build_lib).resolve() / "pygcmc_tools" / "bin"
        package_bin.mkdir(parents=True, exist_ok=True)
        target = package_bin / "gcmc_cpu"
        shutil.copy2(source, target)
        target.chmod(target.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


setup(
    name="pygcmc",
    version="0.1.0",
    description="C++17 Grand Canonical Monte Carlo library with Python bindings",
    long_description=(
        (PROJECT_ROOT / "README.md").read_text(encoding="utf-8")
        if (PROJECT_ROOT / "README.md").exists()
        else ""
    ),
    long_description_content_type="text/markdown",
    ext_modules=[CMakeExtension("pygcmc")],
    cmdclass={"build_ext": CMakeBuild},
    packages=["pygcmc_tools"],
    package_dir={"": "python"},
    package_data={"pygcmc_tools": ["bin/gcmc_cpu"]},
    entry_points={
        "console_scripts": [
            "gcmc_cpu=pygcmc_tools.gcmc_cpu:main",
            "pygcmc_paper_validation=pygcmc_tools.paper_validation:main",
        ]
    },
    python_requires=">=3.9",
    extras_require={
        "test": ["numpy", "pytest", "pytest-timeout", "scipy"],
    },
    zip_safe=False,
)
