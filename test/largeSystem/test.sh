#!/bin/bash

# This script is used to test the functionality of the program.
cd ../../
pip install .
cd test/largeSystem

pygcmc -p memcrowdatp_cg.pdb -s memcrowdatp_cg.psf -o output.txt -e 1024 -n 10

rm charmm36.ff output.txt 
rm -rf ../../pyGCMC.egg-info ../../build ../../gcmc/*.o ../../.eggs ../../.vscode
