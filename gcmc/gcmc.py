"""

    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Abhishek A. Kognole, 
    	Aoxiang Tao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu

"""


# import numpy as np
import pkg_resources
from . import protein_data
import sys
import os

import argparse

# https://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python

class Tee(object):
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush() # If you want the output to be visible immediately
    def flush(self) :
        for f in self.files:
            f.flush()


class GCMC:
    
    def __init__(self):

        # find the path of the mols files        
        self.mols_str_path = pkg_resources.resource_filename(__name__, 'mols') + '/'

        
        self.get_simulation()


        
        
    def get_fragment(self,fragment = None):
        if fragment is None:
            self.fragmentName = ['acey', 'benx', 'dmee', 'form', 'imia', 'mamy', 'meoh', 'prpx', 'sol']


            # fragmentPdb = [pmd.load_file('%s%s.pdb' % (self.mols_str_path,i)) for i in fragmentName]
            # fragmentPsf = [pmd.load_file('%s%s.psf' % (self.mols_str_path,i)) for i in fragmentName]
            
            # for i in range(len(fragmentName)):
            #     fragmentPsf[i].positions = fragmentPdb[i].positions
            # self.fragment = fragmentPsf

            # self.fragment = [(fragmentPdb[i],fragmentPsf[i]) for i in range(len(fragmentName))]


    def get_simulation(self):

        # NBFIX  rmin=<charmm_rmin>/2^(1/6), eps=4.184*<charmm_eps>
        # NB  rmin=<charmm_rmin>/2^(1/6)*2, eps=-4.184*<charmm_eps>

        # self.nbfix_dict = {i:[self.params.nbfix_types[i][1] / 2**(1./6) * 0.1,self.params.nbfix_types[i][0] * 4.184 ] 
        #             for i in self.params.nbfix_types}
        
        # self.nb_dict = {i:[self.params.atom_types_str[i].rmin / 2 ** (1./6) * 0.1 * 2,self.params.atom_types_str[i].epsilon *-4.184] 
        #                 for i in self.params.atom_types_str}
    
        
        
        self.get_fragment()

        

def main():
    # file_output = open('Analyze_output.txt', 'w')
    # original_output = sys.stdout
    # sys.stdout = Tee(sys.stdout, file_output)

    parser = argparse.ArgumentParser(description="pyGCMC - A python package for GCMC simulation")

    parser.add_argument(
        "-P",
        "--pdb-file",
        dest="pdb_file",
        required=True,
        help="The file .pdb for GCMC",
        metavar="file.pdb",
        type=str,
    )
    parser.add_argument(
        "-T",
        "--top-file",
        dest="top_file",
        required=False,
        help="The file .top for GCMC",
        metavar="file.top",
        type=str,
    )
    parser.add_argument(
        "-O",
        "--out-file",
        dest="out_file",
        required=False,
        help="The output file for GCMC",
        metavar="file.txt",
        type=str,
    )
    args = parser.parse_args()
    pdb_file = args.pdb_file
    top_file = args.top_file
    out_file = args.out_file
    if out_file is not None:
        file_output = open(out_file, 'w')
        original_output = sys.stdout
        sys.stdout = Tee(sys.stdout, file_output)
    
    print(f"Using pdb file: {pdb_file}")
    if top_file is not None:
        print(f"Using top file: {top_file}")
        if os.path.exists('temp_link'):
            os.remove('temp_link')
        os.symlink(pkg_resources.resource_filename(__name__, 'charmm36.ff'), 'temp_link') # create a symbolic link to the force field directory
        os.rename('temp_link', 'charmm36.ff') # rename the symbolic link
    print(f"Using output file: {out_file}")

    try:
        cryst, atoms = protein_data.read_pdb(pdb_file)
    except:
        print(f"Error reading pdb file: {pdb_file}")
        sys.exit(1)
    print(f"pdb cryst: {cryst}")
    print(f"pdb atom number: {len(atoms)}")

    protein_data.read_top(top_file)

    # if top_file is not None:
    #     try:
    #         # top = protein_data.read_top(top_file)
    #         protein_data.read_top(top_file)
    #     except:
    #         print(f"Error reading top file: {top_file}")
    #         sys.exit(1)
    #     # print(f"top atom number: {len(top)}")

    if out_file is not None:
        sys.stdout = original_output
        file_output.close()