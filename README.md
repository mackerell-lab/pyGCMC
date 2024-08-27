# pyGCMC

pyGCMC is a python package for performing grand canonical Monte Carlo simulations of molecules in nanoporous materials. It is designed to be easy to use and to be easily extensible. The package is written in Python, C++, and CUDA.

Current version: 1.1.240324

## Requirements

- Python 3.6 or later
- numpy 1.18 or later
- C++11 compatible compiler
- CUDA compatible GPU and CUDA toolkit (version X.X or later)

## Installation

The package can be downloaded from GitHub https://github.com/mackerell-lab/pyGCMC

Download the source code from GitHub and run the setup script:
```
pip install .
```
Or you can install the package from PyPI:
```
pip install pygcmc
```
The web page for the package on PyPI is https://pypi.org/project/pygcmc/

## Usage

pyGCMC provides two command-line interfaces: `pygcmc` (new version) and `gcmc` (old version).

### New Version (pygcmc)

Run the pyGCMC executable:
```
usage: pygcmc [-h] -p file.pdb [-t file.top] [-s file.psf] [-o file.txt] [-u muex1,muex2,...] [-f conf1,conf2,... or conf] [-n mcsteps] [-m mctime1,mctime2,...]
            [-c conc1,conc2,...] [-y cavity_bias_dx] [-e seed] [-P] [-w]

pyGCMC - A python package for GCMC simulation

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
                        The value of solute muex(splice by , with no space), if the first value is negative, then follow the -u or --fragmuex without space
  -f conf1,conf2,... or conf, --fragconf conf1,conf2,... or conf
                        The value of solute conf(splice by , with no space). Or only one value for all solutes
  -n mcsteps, --mcsteps mcsteps
                        The number of MC steps
  -m mctime1,mctime2,..., --mctime mctime1,mctime2,...
                        The mctime of solutes(splice by , with no space)
  -c conc1,conc2,..., --fragconc conc1,conc2,...
                        The value of solute concentration(splice by , with no space)
  -y cavity_bias_dx, --cavitybias-dx cavity_bias_dx
                        The value of cavity bias dx(if dx <= 0, then no cavity bias)
  -e seed, --seed seed  The seed of random number
  -P, --PME             Enable PME(Default: Disable)
  -w, --show-info       Show the information of solutes

  ```

### Old Version (gcmc)

To run `gcmc` as the old version, you can use the following command:

```
usage: gcmc [-h] -p PARAMFILE [-v] [--logfile LOGFILE] [--debug] [--version]

Perform GCMC Simulation

options:
  -h, --help            show this help message and exit
  -p PARAMFILE, --paramfile PARAMFILE
                        [Required] input parameter file
  -v, --verbose         [Optional] verbose output
  --logfile LOGFILE     [Optional] log file, if not specified, then output will be stdout
  --debug               [Optional] for debug purpose
  --version             Show program's version number and exit
  ```
