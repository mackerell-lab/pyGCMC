#!/bin/bash

# This script is used to test the functionality of the program.
cd ..
pip install .
cd test
#for ((i=0;i<100;i++))
#do
#echo "seed = $i" &>> tmp~
gcmc -p 6v3g_silcs.1.pdb -t 6v3g_silcs.1.top -o output.txt -e 10 
#done

rm charmm36.ff output.txt 
rm -rf ../gcmc.egg-info ../build ../gcmc/*.o ../.eggs ../.vscode
