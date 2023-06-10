"""

    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Abhishek A. Kognole, 
    	Aoxiang Tao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu

"""


import numpy as np
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

        self.fragmentName = ['BENX', 'PRPX', 'DMEE', 'MEOH', 'FORM', 'IMIA', 'ACEY', 'MAMY', 'SOL']

        self.ff_files = ['charmm36.ff/ffnonbonded.itp', 'charmm36.ff/nbfix.itp', 'charmm36.ff/silcs.itp']

        self.top_file = None
        self.pdb_file = None

        if os.path.exists('temp_link'):
            os.remove('temp_link')
        os.symlink(pkg_resources.resource_filename(__name__, 'charmm36.ff'), 'temp_link') # create a symbolic link to the force field directory
        os.rename('temp_link', 'charmm36.ff') # rename the symbolic link to force field directory




    
    def get_pdb(self,pdb_file):
        self.pdb_file = pdb_file

        print(f"Using pdb file: {pdb_file}")


        try:
            self.cryst, self.atoms = protein_data.read_pdb(pdb_file)
        except:
            print(f"Error reading pdb file: {pdb_file}")
            sys.exit(1)


        print(f"pdb cryst: {self.cryst}")
        print(f"pdb atom number: {len(self.atoms)}")
            
    def get_top(self,top_file):
        self.top_file = top_file


        print(f"Using top file: {top_file}")
        

        try:
            self.atom_top = protein_data.read_top(top_file)
        except:
            print(f"Error reading top file: {top_file}")
            sys.exit(1)

        print(f"top atom number: {len(self.atom_top)}")

        if len(self.atoms) != len(self.atom_top):
            print(f"Error: pdb atom number {len(self.atoms)} != top atom number {len(self.atom_top)}")
            sys.exit(1)
    
        for i, atom in enumerate(self.atoms):
            # if atom.name != atom_top[i][3]:
            #     print(f"Error: pdb atom {i+1} name {atom.name} != top atom name {atom_top[i][3]}")
            #     sys.exit(1)
            if atom.residue[:3] != self.atom_top[i][2][:3]:
                print(f"Error: pdb atom {i+1} residue {atom.residue} != top atom residue {self.atom_top[i][2]}")
                sys.exit(1)
            atom.residue = self.atom_top[i][2]
            atom.type = self.atom_top[i][0]
            atom.charge = self.atom_top[i][4]
            atom.nameTop = self.atom_top[i][3]
        


    def get_fragment(self):
        self.fragments = []

        for frag in self.fragmentName:

            try:

                self.fragments += [protein_data.read_pdb(f"charmm36.ff/mol/{frag.lower()}.pdb")[1]]
                atom_top = protein_data.read_itp(f"charmm36.ff/mol/{frag.lower()}.itp")
                for i, atom in enumerate(self.fragments[-1]):
                    # if atom.name != atom_top[i][3]:
                    #     print(f"Error: pdb atom {i+1} name {atom.name} != top atom name {atom_top[i][3]}")
                    #     sys.exit(1)
                    if atom.residue[:3] != atom_top[i][2][:3]:
                        print(f"Error: pdb atom {i+1} residue {atom.residue} != top atom residue {atom_top[i][2]}")
                        sys.exit(1)
                    atom.residue = atom_top[i][2]
                    atom.type = atom_top[i][0]
                    atom.charge = atom_top[i][4]
                    atom.nameTop = atom_top[i][3]
            except:
                print(f"Error reading fragment file: charmm36.ff/mol/{frag}")
                sys.exit(1)

            print(f"Read fragment: {frag} with {len(self.fragments[-1])} atoms")
    # def get_fragment(self,fragment = None):
    #     if fragment is None:
    #         self.fragmentName = ['BENX', 'PRPX', 'DMEE', 'MEOH', 'FORM', 'IMIA', 'ACEY', 'MAMY', 'SOL']


    def get_forcefield(self):
        self.nb_dict = {}
        self.nbfix_dict = {}

        for ff_file in self.ff_files:
            try:
                nb_dict, nbfix_dict = protein_data.read_ff(ff_file)
                print(f"Using force field file: {ff_file} ,\tnb number: {len(nb_dict)} ,\tnbfix number: {len(nbfix_dict)}")
                self.nb_dict.update(nb_dict)
                self.nbfix_dict.update(nbfix_dict)
            except:
                print(f"Error reading force field file: {ff_file}")
                sys.exit(1)
    
    

    def get_simulation(self):

        
        self.get_fragment()
        self.get_forcefield()
        
        if self.pdb_file is None or self.top_file is None:
            print("Error: pdb or top file not set")
            sys.exit(1)
        
        self.atomtypes1 = []
        self.atomtypes2 = []
        for frag in self.fragments:
            for atom in frag:
                if atom.type not in self.atomtypes2:
                    self.atomtypes1.append(atom.type)
                    self.atomtypes2.append(atom.type)
        
        for atom in self.atoms:
            if atom.type not in self.atomtypes2:
                self.atomtypes2.append(atom.type)
        
        print(f"Fragment Type number: {len(self.atomtypes1)}")
        print(f"Total Type number: {len(self.atomtypes2)}")

        self.ff_pairs = []

        for i, type1 in enumerate(self.atomtypes1):
            for j, type2 in enumerate(self.atomtypes2):
                ff_pair = [0,0]
                if (type1, type2) in self.nbfix_dict:
                    ff_pair[0] = self.nbfix_dict[(type1, type2)][0]
                    ff_pair[1] = self.nbfix_dict[(type1, type2)][1]
                elif (type2, type1) in self.nbfix_dict:
                    ff_pair[0] = self.nbfix_dict[(type2, type1)][0]
                    ff_pair[1] = self.nbfix_dict[(type2, type1)][1]
                else:
                    sigma1 = self.nb_dict[type1][0]
                    sigma2 = self.nb_dict[type2][0]
                    sigma = (sigma1 + sigma2) / 2
                    epsilon1 = self.nb_dict[type1][1]
                    epsilon2 = self.nb_dict[type2][1]
                    epsilon = np.sqrt(epsilon1 * epsilon2)

                    ff_pair[0] = sigma
                    ff_pair[1] = epsilon
                self.ff_pairs.append(ff_pair)

        
        print(f"Total FF pair number: {len(self.ff_pairs)}")
        
        # print(self.ff_pairs)

        for atom in self.atoms:
            atom.typeNum = self.atomtypes2.index(atom.type)
            
        for frag in self.fragments:
            for atom in frag:
                atom.typeNum = self.atomtypes2.index(atom.type)


        self.fix_atoms = []
        relax_atoms = [[] for i in range(len(self.fragments))]
        for atom in self.atoms:
            if atom.residue in self.fragmentName:
                relax_atoms[self.fragmentName.index(atom.residue)].append(atom)
            else:
                self.fix_atoms.append(atom)
        

        self.fraglist =[[relax_atoms[i][j:j + len(self.fragments[i])] for j in range(0, len(relax_atoms[i]), len(self.fragments[i]))] for i in range(len(self.fragments))]
        
        print(f"Total fixed atom number: {len(self.fix_atoms)}")
        print(f"Total relax atom number: {len(self.atoms) - len(self.fix_atoms)}")
        
        for i, frag in enumerate(self.fragments):
            print(f"Fragment {self.fragmentName[i]} number: {len(self.fraglist[i])}")

        
        # for atom in self.fix_atoms:
        #     print(f"Fix atom: {atom.name} {atom.nameTop} {atom.type} {atom.typeNum}")


        
        # for frag in self.fragments:
        #     print(f"Fragment: {frag[0].residue}")
        #     for atom in frag:
        #         print(f"{atom.name} {atom.type} {atom.typeNum}")

        # self.residueType = []
        # for frag in self.fragments:
        #     for atom in frag:
        #         if atom.residue not in self.residueType:
        #             self.residueType.append(atom.residue)
        
        # for atom in self.atoms:
        #     if atom.residue not in self.residueType:
        #         self.residueType.append(atom.residue)

        # for i, residue in enumerate(self.residueType):
        #     print(f"Residue {i}: {residue}")




    # def get_simulation(self):

    #     # NBFIX  rmin=<charmm_rmin>/2^(1/6), eps=4.184*<charmm_eps>
    #     # NB  rmin=<charmm_rmin>/2^(1/6)*2, eps=-4.184*<charmm_eps>

    #     # self.nbfix_dict = {i:[self.params.nbfix_types[i][1] / 2**(1./6) * 0.1,self.params.nbfix_types[i][0] * 4.184 ] 
    #     #             for i in self.params.nbfix_types}
        
    #     # self.nb_dict = {i:[self.params.atom_types_str[i].rmin / 2 ** (1./6) * 0.1 * 2,self.params.atom_types_str[i].epsilon *-4.184] 
    #     #                 for i in self.params.atom_types_str}
    
        
        
    

        

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
    
    gcmc = GCMC()

    pdb_file = args.pdb_file
    top_file = args.top_file
    out_file = args.out_file
    if out_file is not None:
        file_output = open(out_file, 'w')
        original_output = sys.stdout
        sys.stdout = Tee(sys.stdout, file_output)
        print(f"Using output file: {out_file}")
    

    gcmc.get_pdb(pdb_file)

        
    




    if top_file is not None:
        gcmc.get_top(top_file)

    
    gcmc.get_simulation()

    


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