from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import os
import shutil
import subprocess
# import numpy as np

# CUDA specific config

nvcc_bins = [os.environ.get('OPENMM_CUDA_COMPILER'),shutil.which('nvcc'),'/usr/local/cuda/bin/nvcc']

nvcc_bin = [nvcc_path for nvcc_path in nvcc_bins if nvcc_path and os.path.exists(nvcc_path)][0]

nvcc_dir = os.path.dirname(os.path.abspath(nvcc_bin))

cuda_lib_path = os.path.join(os.path.dirname(nvcc_dir), 'lib64')

cuda_include_path = os.path.join(os.path.dirname(nvcc_dir), 'include')

print("\nnvcc_bin:\t\t", nvcc_bin,"\ncuda_lib_path:\t\t", cuda_lib_path,"\ncuda_include_path:\t", cuda_include_path,"\n")

package_data = {
    'gcmc': ['toppar.str', 'resources.zip', 'charmm36.ff/*','toppar/*','charmm36.ff/mol/*','*.cu','*.h','*.cpp'],
}

# Load README.md for the long description

with open("README.md", "r") as fh:
    long_description = fh.read()


class CustomBuildExt(build_ext):
    def build_extensions(self):
        import numpy as np

        _gcc_compile_args = ["-std=c++11", "-fPIC"]
        
        # Compile the CUDA code
        cuda_file = "gcmc/gcmc.cu"
        obj_file = "gcmc/gcmc.o"
        nvcc_command = [nvcc_bin, "-c", cuda_file, "-o", obj_file, "--compiler-options", "-fPIC"]
        subprocess.check_call(nvcc_command)

        for ext in self.extensions:
            ext.extra_compile_args = _gcc_compile_args
            ext.extra_objects = [obj_file]  # Link the CUDA object file
            ext.include_dirs.append(cuda_include_path)  # Add CUDA include path
            ext.include_dirs.append(np.get_include())
            # ext.include_dirs.append(os.path.join(os.path.dirname(__file__), 'gcmc')) 
            ext.library_dirs.append(cuda_lib_path)  # Add CUDA library path
            ext.libraries.append("cudart")  # Add the CUDA runtime library
        super().build_extensions()

ext_modules = [
    Extension(
        "gcmc.gpu",
        sources=["gcmc/gcmc.cpp"],
        language="c++",
    )
]

setup(
    install_requires=["numpy"],
    setup_requires=["numpy"],
    name="pyGCMC",
    version="0.9.231216",
    packages=find_packages(),
    package_data=package_data,
    ext_modules=ext_modules,
    cmdclass={"build_ext": CustomBuildExt},
    entry_points={
        'console_scripts': [
            'pygcmc=gcmc:main',
            'gcmc=gcmc:mainOld'
        ],
    },
    long_description=long_description,
    long_description_content_type='text/markdown'
)
