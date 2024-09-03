"""
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu
"""

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from distutils.command.clean import clean
import os
import shutil
import subprocess
import platform

def find_nvcc():
    nvcc_bins = [
        os.environ.get('OPENMM_CUDA_COMPILER'),
        shutil.which('nvcc'),
        '/usr/local/cuda/bin/nvcc'
    ]
    for nvcc_path in nvcc_bins:
        if nvcc_path and os.path.exists(nvcc_path):
            return nvcc_path
    raise RuntimeError("CUDA compiler (nvcc) not found. Please ensure CUDA is installed and in your PATH.")

def get_cuda_version(nvcc_path):
    try:
        result = subprocess.run([nvcc_path, '--version'], capture_output=True, text=True)
        version_str = result.stdout.split('release ')[-1].split(',')[0]
        return tuple(map(int, version_str.split('.')))
    except:
        return None

# CUDA specific configuration
nvcc_bin = find_nvcc()
cuda_version = get_cuda_version(nvcc_bin)

if cuda_version and cuda_version < (9, 0):
    raise RuntimeError(f"CUDA version {'.'.join(map(str, cuda_version))} is not supported. Please use CUDA 9.0 or later.")

nvcc_dir = os.path.dirname(os.path.abspath(nvcc_bin))
if platform.system() == 'Windows':
    cuda_lib_path = os.path.join(os.path.dirname(nvcc_dir), 'lib', 'x64')
else:
    cuda_lib_path = os.path.join(os.path.dirname(nvcc_dir), 'lib64')
cuda_include_path = os.path.join(os.path.dirname(nvcc_dir), 'include')

print(f"\nnvcc_bin: {nvcc_bin}")
print(f"cuda_lib_path: {cuda_lib_path}")
print(f"cuda_include_path: {cuda_include_path}")
print(f"CUDA Version: {'.'.join(map(str, cuda_version)) if cuda_version else 'Unknown'}\n")

# Package data
package_data = {
    'gcmc': [
        'toppar.str', 
        'resources.zip', 
        'charmm36.ff/*', 
        'toppar/*', 
        'charmm36.ff/mol/*', 
        '*.cu', 
        '*.h', 
        '*.cpp'
    ],
}

# Load README.md for the long description
with open("README.md", "r", encoding="utf-8") as file:
    long_description = file.read()

class CustomBuildExt(build_ext):
    def build_extensions(self):
        import numpy as np

        gcc_compile_args = ["-std=c++11", "-fPIC", "-O3"]
        
        # Compile the CUDA code
        cuda_file = "gcmc/gcmc.cu"
        obj_file = "gcmc/gcmc.o"
        nvcc_command = [nvcc_bin, "-c", cuda_file, "-o", obj_file, "--compiler-options", "-fPIC", "-O3"]
        subprocess.check_call(nvcc_command)

        for ext in self.extensions:
            ext.extra_compile_args = gcc_compile_args
            ext.extra_objects = [obj_file]
            ext.include_dirs.extend([cuda_include_path, np.get_include()])
            ext.library_dirs.append(cuda_lib_path)
            ext.libraries.append("cudart")
        super().build_extensions()

class CustomClean(clean):
    def run(self):
        super().run()
        files_to_delete = ['gcmc/gcmc.o']
        for file in files_to_delete:
            try:
                os.remove(file)
                print(f"Removed {file}")
            except OSError:
                pass

# Extension modules
ext_modules = [
    Extension(
        "gcmc.gpu",
        sources=["gcmc/gcmc.cpp"],
        language="c++",
    )
]

# Setup configuration
setup(
    name="pyGCMC",
    version="1.3.240903",
    author="Mingtian Zhao, Alexander D. MacKerell Jr.",
    author_email="zhaomt@outerbanks.umaryland.edu",
    description="A python package for performing grand canonical Monte Carlo simulations",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/mackerell-lab/pyGCMC",
    packages=find_packages(),
    package_data=package_data,
    ext_modules=ext_modules,
    cmdclass={"build_ext": CustomBuildExt, "clean": CustomClean},
    entry_points={
        'console_scripts': [
            'pygcmc=gcmc:main',
            'gcmc=gcmc:mainOld'
        ],
    },
    install_requires=["numpy>=1.18,<2"],
    setup_requires=["numpy>=1.18,<2"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)