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
import random

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

        self.ff_files = ['charmm36.ff/ffnonbonded.itp', 'charmm36.ff/nbfix.itp', 'charmm36.ff/silcs.itp']


        self.fragmentName = ['BENX', 'PRPX', 'DMEE', 'MEOH', 'FORM', 'IMIA', 'ACEY', 'MAMY', 'SOL']

        self.fragconc = [ 0.25,   0.25,   0.25,   0.25,   0.25,   0.25,   0.25,   0.25,  55.00]

        self.fragmuex = [-2.79, 1.46, -1.44, -5.36, -10.92, -14.18, -97.31, -68.49, -5.6]

        self.fragconf = [1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]

        self.mctime = [1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000]

        self.mcsteps = 10000


        self.attempt_prob_frag = [0.300, 0.300, 0.200, 0.200]
        
        self.attempt_prob_water = [0.250, 0.250, 0.250, 0.250]



        self.cutoff = 15.0

        self.fixCutoff = 6.0
        
        self.grid_dx = 1.0


        
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

        print(f"pdb atom number: {len(self.atoms)}")   

        print(f"pdb cryst: {self.cryst}")
        

        x_values = [atom.x for atom in self.atoms]
        y_values = [atom.y for atom in self.atoms]
        z_values = [atom.z for atom in self.atoms]

        

        print("x direction: min =", min(x_values), "max =", max(x_values))
        print("y direction: min =", min(y_values), "max =", max(y_values))
        print("z direction: min =", min(z_values), "max =", max(z_values))

        
            
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
    
    def get_move(self):
        self.move_array = random.choices(range(len(self.fragmentName)), weights=self.mctime, k=self.mcsteps)

        for i, n in enumerate(self.move_array):
            if self.fragmentName[n] == 'SOL':
                self.move_array[i] = self.move_array[i] * 4 + random.choices(range(len(self.attempt_prob_water)), weights=self.attempt_prob_water, k=1)[0]
            else:
                self.move_array[i] = self.move_array[i] * 4 + random.choices(range(len(self.attempt_prob_frag)), weights=self.attempt_prob_frag, k=1)[0]

        self.move_array_n = [0 for i in range(len(self.attempt_prob_frag) * len(self.fragmentName))]
        self.move_array_frag = [0 for i in range(len(self.fragmentName))]

        for i in self.move_array:
            self.move_array_n[i] += 1
            self.move_array_frag[i//4] += 1
        
        for i in range(len(self.fragmentName)):
            print(f"Fragment {self.fragmentName[i]}\t move {self.move_array_frag[i]} times", end='\t')
            for j in range(len(self.attempt_prob_frag)):
                print(f"Movement {j} move {self.move_array_n[i*4+j]} times", end='\t')
            print()

        # for i in range(len(self.attempt_prob_frag)):
        #     print(f"Movement {i} move {random_n[i]} times")
        
        # random_n = [0 for i in range(len(self.attempt_prob_frag))]

        # for i in random_array:
        #     random_n[i%4] += 1
        
        # for i, frag in enumerate(random_n):
        #     print(f"Movement {i} move {random_n[i]} times")


    # def update_data(self):
    #     pass

    def get_fragmuex(self, fragmuex):
        
        fragmuex = [float(i) for i in fragmuex if len(i) > 0]
        if len(fragmuex) != len(self.fragmentName):
            print("Error: fragmuex number not match")
            sys.exit(1)
        else:
            self.fragmuex = fragmuex
    
    def get_fragconf(self, fragconf):

        fragconf = [int(i) for i in fragconf if len(i) > 0]
        if len(fragconf) == 1:
            fragconf = fragconf * len(self.fragmentName)
            self.fragconf = fragconf
        elif len(fragconf) == len(self.fragmentName):
            self.fragconf = fragconf
        else:
            print("Error: fragconf number not match")
            sys.exit(1)

    def get_mctime(self, mctime):

        mctime = [float(i) for i in mctime if len(i) > 0]
        if len(mctime) != len(self.fragmentName):
            print("Error: mctime number not match")
            sys.exit(1)
        else:
            self.mctime = mctime

    def get_fragconc(self, fragconc):

        fragconc = [float(i) for i in fragconc if len(i) > 0]
        if len(fragconc) != len(self.fragmentName):
            print("Error: fragconc number not match")
            sys.exit(1)
        else:
            self.fragconc = fragconc


    def show_parameters(self):
        print(f"MC steps: {self.mcsteps}")
        print("Fragment Name: \t\t",'\t\t'.join(self.fragmentName))
        print("Fragment Muex: \t\t",'\t\t'.join([str(i) for i in self.fragmuex]))
        print("Fragment Conc: \t\t",'\t\t'.join([str(i) for i in self.fragconc]))
        print("Fragment ConfBias: \t",'\t\t'.join([str(i) for i in self.fragconf]))
        print("Fragment mcTime: \t",'\t\t'.join([str(i) for i in self.mctime]))


        


    def get_simulation(self):

        
        self.get_fragment()
        self.get_forcefield()

        self.get_move()

        self.show_parameters()
        
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


        self.set_fixCut()

        # self.update_data()

        # for atom in self.fix_atoms:
        #     print(f"Fix atom: {atom.name} {atom.nameTop} {atom.type} {atom.typeNum} {atom.sequence} {atom.sequence2}")

        

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

        # atoms_list = []
        # atoms_sec = None
        # for atom in self.fix_atoms:
        #     if atom.sequence != atoms_sec:
        #         atoms_list.append([atom])
        #         atoms_sec = atom.sequence
        #     else:
        #         atoms_list[-1].append(atom)
        # print(f"Total residue number: {len(atoms_list)}")


        # def distance(atom1, atom2):
        #     x1, y1, z1 = atom1.x, atom1.y, atom1.z
        #     x2, y2, z2 = atom2.x, atom2.y, atom2.z
            
        #     return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

        # def find_max_distance(atom_list):
        #     max_dist = 0
            
        #     for i, atom1 in enumerate(atom_list[:-1]):
        #         for atom2 in atom_list[i + 1:]:
        #             dist = distance(atom1, atom2)
        #             max_dist = max(max_dist, dist)
            
        #     return max_dist

        # for atom_list in atoms_list:
        #     max_dist = find_max_distance(atom_list)
        #     print(f"Max distance: {max_dist:.3f}")




    # def get_simulation(self):

    #     # NBFIX  rmin=<charmm_rmin>/2^(1/6), eps=4.184*<charmm_eps>
    #     # NB  rmin=<charmm_rmin>/2^(1/6)*2, eps=-4.184*<charmm_eps>

    #     # self.nbfix_dict = {i:[self.params.nbfix_types[i][1] / 2**(1./6) * 0.1,self.params.nbfix_types[i][0] * 4.184 ] 
    #     #             for i in self.params.nbfix_types}
        
    #     # self.nb_dict = {i:[self.params.atom_types_str[i].rmin / 2 ** (1./6) * 0.1 * 2,self.params.atom_types_str[i].epsilon *-4.184] 
    #     #                 for i in self.params.atom_types_str}
    
        
        
    
    def set_fixCut(self):
        
        x0 = float('inf')
        y0 = float('inf')
        z0 = float('inf')
        n = -1
        for atom in self.fix_atoms:
            x, y, z = atom.x, atom.y, atom.z
            if np.sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2) > self.fixCutoff:
                n += 1
                x0, y0, z0 = x, y, z
                atom.sequence2 = n
            else:
                atom.sequence2 = n
        

def main():
    # file_output = open('Analyze_output.txt', 'w')
    # original_output = sys.stdout
    # sys.stdout = Tee(sys.stdout, file_output)

    parser = argparse.ArgumentParser(description="pyGCMC - A python package for GCMC simulation")

    parser.add_argument(
        "-p",
        "--pdb-file",
        dest="pdb_file",
        required=True,
        help="The file .pdb for GCMC",
        metavar="file.pdb",
        type=str,
    )
    parser.add_argument(
        "-t",
        "--top-file",
        dest="top_file",
        required=False,
        help="The file .top for GCMC",
        metavar="file.top",
        type=str,
    )
    parser.add_argument(
        "-o",
        "--out-file",
        dest="out_file",
        required=False,
        help="The output file for GCMC",
        metavar="file.txt",
        type=str,
    )
    parser.add_argument(
        "-u",
        "--fragmuex",
        dest="fragmuex",
        required=False,
        help="The value of fragment muex(splice by , with no space)",
        metavar="muex1,muex2,...",
        type=str,
    ) 
    parser.add_argument(
        "-f",
        "--fragconf",
        dest="fragconf",
        required=False,
        help="The value of fragment conf(splice by , with no space). Or only one value for all fragments",
        metavar="conf1,conf2,... or conf",
        type=str,
    )
    parser.add_argument(
        "-n",
        "--mcsteps",
        dest="mcsteps",
        required=False,
        help="The number of MC steps",
        metavar="mcsteps",
        type=int,
    )
    parser.add_argument(
        "-m",
        "--mctime",
        dest="mctime",
        required=False,
        help="The mctime of Fragments(splice by , with no space)",
        metavar="mctime1,mctime2,...",
        type=str,
    )
    parser.add_argument(
        "-c",
        "--fragconc",
        dest="fragconc",
        required=False,
        help="The value of fragment concentration(splice by , with no space)",
        metavar="conc1,conc2,...",
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

    if args.fragmuex is not None:
        fragmuex = args.fragmuex
        fragmuex = fragmuex.split(',')
        gcmc.get_fragmuex(fragmuex)
        print(f"Using fragment muex: {fragmuex}")    

    if args.fragconf is not None:
        fragconf = args.fragconf
        fragconf = fragconf.split(',')
        gcmc.get_fragconf(fragconf)
        print(f"Using fragment conf: {fragconf}")

    if args.mcsteps is not None:
        gcmc.mcsteps = args.mcsteps
        print(f"Using MC steps: {args.mcsteps}")

    if args.mctime is not None:
        mctime = args.mctime
        mctime = mctime.split(',')
        gcmc.get_mctime(mctime)
        print(f"Using fragment mctime: {mctime}")

    if args.fragconc is not None:
        fragconc = args.fragconc
        fragconc = fragconc.split(',')
        gcmc.get_fragconc(fragconc)
        print(f"Using fragment concentration: {fragconc}")

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