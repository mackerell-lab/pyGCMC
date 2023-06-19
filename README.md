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
usage: gcmc [-h] -p file.pdb [-t file.top] [-s file.psf] [-o file.txt] [-u muex1,muex2,...] [-f conf1,conf2,... or conf] [-n mcsteps] [-m mctime1,mctime2,...] [-c conc1,conc2,...] [-y cavity_bias_dx] [-e seed] [-w]

options:
  -h, --help            show this help message and exit
  -p file.pdb, --pdb-file file.pdb
                        The file .pdb for GCMC
  -t file.top, --top-file file.top
                        The file .top for GCMC
  -s file.psf, --psf-file file.psf
                        The file .psf for GCMC
  -o file.txt, --out-file file.txt
                        The output file for GCMC
  -u muex1,muex2,..., --fragmuex muex1,muex2,...
                        The value of fragment muex(splice by , with no space)
  -f conf1,conf2,... or conf, --fragconf conf1,conf2,... or conf
                        The value of fragment conf(splice by , with no space). Or only one value for all fragments
  -n mcsteps, --mcsteps mcsteps
                        The number of MC steps
  -m mctime1,mctime2,..., --mctime mctime1,mctime2,...
                        The mctime of Fragments(splice by , with no space)
  -c conc1,conc2,..., --fragconc conc1,conc2,...
                        The value of fragment concentration(splice by , with no space)
  -y cavity_bias_dx, --cavitybias-dx cavity_bias_dx
                        The value of cavity bias dx(if dx <= 0, then no cavity bias)
  -e seed, --seed seed  The seed of random number
  -w, --show-info       Show the information of fragments

  ```

