"""

    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu

"""

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
import os
import shutil
import subprocess

# CUDA specific configuration
nvccBins = [os.environ.get('OPENMM_CUDA_COMPILER'), shutil.which('nvcc'), '/usr/local/cuda/bin/nvcc']
nvccBin = next((nvccPath for nvccPath in nvccBins if nvccPath and os.path.exists(nvccPath)), None)

nvccDir = os.path.dirname(os.path.abspath(nvccBin))
cudaLibPath = os.path.join(os.path.dirname(nvccDir), 'lib64')
cudaIncludePath = os.path.join(os.path.dirname(nvccDir), 'include')

print("\nnvcc_bin:\t\t", nvccBin, 
      "\ncuda_lib_path:\t\t", cudaLibPath, 
      "\ncuda_include_path:\t", cudaIncludePath, "\n")

# Package data
packageData = {
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
with open("README.md", "r") as file:
    longDescription = file.read()


# Custom build extension class
class CustomBuildExt(build_ext):
    """ Custom build extension class. """
    def build_extensions(self):
        import numpy as np

        gccCompileArgs = ["-std=c++11", "-fPIC"]
        
        # Compile the CUDA code
        cudaFile = "gcmc/gcmc.cu"
        objFile = "gcmc/gcmc.o"
        nvccCommand = [nvccBin, "-c", cudaFile, "-o", objFile, "--compiler-options", "-fPIC"]
        subprocess.check_call(nvccCommand)

        for ext in self.extensions:
            ext.extra_compile_args = gccCompileArgs
            ext.extra_objects = [objFile]  # Link the CUDA object file
            ext.include_dirs.append(cudaIncludePath)  # Add CUDA include path
            ext.include_dirs.append(np.get_include())
            ext.library_dirs.append(cudaLibPath)  # Add CUDA library path
            ext.libraries.append("cudart")  # Add the CUDA runtime library
        super().build_extensions()

# Extension modules
extModules = [
    Extension(
        "gcmc.gpu",
        sources=["gcmc/gcmc.cpp"],
        language="c++",
    )
]

# Setup configuration
setup(
    install_requires=["numpy"],
    setup_requires=["numpy"],
    name="pyGCMC",
    version="1.1.240324",
    packages=find_packages(),
    package_data=packageData,
    ext_modules=extModules,
    cmdclass={"build_ext": CustomBuildExt},
    entry_points={
        'console_scripts': [
            'pygcmc=gcmc:main',
            'gcmc=gcmc:mainOld'
        ],
    },
    long_description=longDescription,
    long_description_content_type='text/markdown'
)
