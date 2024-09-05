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
import glob

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
        'resources/toppar.str', 
        'resources/resources.zip', 
        'resources/charmm36.ff/*', 
        'resources/toppar/*', 
        'resources/charmm36.ff/mol/*', 
        'cpp/*.cu', 
        'cpp/*.h', 
        'cpp/*.cpp',
        'scripts/*.py',
        'python/*.py'
    ],
}

# Load README.md for the long description
with open("README.md", "r", encoding="utf-8") as file:
    long_description = file.read()

class CustomBuildExt(build_ext):
    def build_extensions(self):
        import numpy as np

        compiler = self.compiler.compiler_type

        if compiler == 'msvc':  # Microsoft Visual C++
            compile_args = ['/O2']
            if platform.machine().endswith('64'):
                compile_args.append('/arch:AVX2')
        else:  # GCC or Clang
            compile_args = ['-O3', '-fPIC', '-march=native']
            if compiler == 'unix':
                compile_args.append('-std=c++11')

        # Compile the CUDA code
        cuda_file = "gcmc/cpp/gcmc.cu"
        obj_file = "gcmc/cpp/gcmc.o"
        nvcc_command = [nvcc_bin, "-c", cuda_file, "-o", obj_file, "--compiler-options", "-fPIC", "-O3"]
        if platform.system() != 'Windows':
            nvcc_command.extend(["-ccbin", "g++"])
        subprocess.check_call(nvcc_command)

        for ext in self.extensions:
            ext.extra_compile_args = compile_args
            ext.extra_objects = [obj_file]
            ext.include_dirs.extend([cuda_include_path, np.get_include()])
            ext.library_dirs.append(cuda_lib_path)
            ext.libraries.append("cudart")
        super().build_extensions()

class CustomClean(clean):
    def run(self):
        super().run()
        files_to_delete = ['gcmc/cpp/gcmc.o', 'gcmc/cpp/*.so', 'gcmc/cpp/*.pyd', 'build', 'pyGCMC.egg-info','.eggs']
        for file_pattern in files_to_delete:
            for file in glob.glob(file_pattern):
                try:
                    if os.path.isfile(file):
                        os.remove(file)
                    elif os.path.isdir(file):
                        shutil.rmtree(file)
                except OSError:
                    pass

# Extension modules
ext_modules = [
    Extension(
        "gcmc.gpu",
        sources=["gcmc/cpp/gcmc.cpp"],
        language="c++",
    )
]

# Setup configuration
setup(
    name="pyGCMC",
    version="1.3.240905",
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
            'pygcmc=gcmc.scripts.main:main',
            'gcmc=gcmc.scripts.main_old:main'
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