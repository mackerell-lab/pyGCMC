# pyGCMC
pyGCMC is a python package for performing grand canonical Monte Carlo simulations of molecules in nanoporous materials. It is designed to be easy to use and to be easily extensible. It is also designed to be easy to install, with minimal dependencies. The package is written in pure python, and the only dependencies are numpy.

## Requirements
The python package is written in python, C++ and CUDA. The python package requires numpy. The C++ code requires a C++11 compatible compiler. The CUDA code requires a CUDA compatible GPU and the CUDA toolkit.

## Installation
Download the source code from GitHub and run the setup script:
```
pip install .
```

## Usage
Run the pyGCMC executable:
```
gcmc [-h] -p file.pdb [-t file.top] [-o file.txt] [-m muex1,muex2,...] [-c conf1,conf2,...] [-n mcsteps]
```

