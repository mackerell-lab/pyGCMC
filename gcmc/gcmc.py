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

import time

from .gpu import runGCMC


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

MOLES_TO_MOLECULES = 0.0006023

Atom_dtype = np.dtype([
    ('position', np.float32, (3, )),
    ('charge', np.float32),
    ('type', np.int32)
])

AtomArray_dtype = np.dtype([
    ('name', 'S4'),

    ('startRes', np.int32),

    ('muex', np.float32),
    ('conc', np.float32),
    ('confBias', np.int32),
    ('mcTime', np.float32),

    ('totalNum', np.int32),
    ('maxNum', np.int32),

    ('num_atoms', np.int32),
    ('atoms', Atom_dtype, (20, ))
])

Info_dtype = np.dtype([
    ('mcsteps', np.int32),
    ('cutoff', np.float32),
    ('grid_dx', np.float32),
    ('startxyz', np.float32, (3, )),
    ('cryst', np.float32, (3, )),
    ('showInfo', np.int32),
    ('cavityFactor', np.float32),
    ('fragTypeNum', np.int32),
    ('totalGridNum', np.int32),
    ('totalResNum', np.int32),
    ('totalAtomNum', np.int32),
    ('ffXNum', np.int32),
    ('ffYNum', np.int32),
    ('PME', np.int32),
    ('seed', np.uint32)
])

Residue_dtype = np.dtype([
    ('position', np.float32, (3, )),
    ('atomNum', np.int32),
    ('atomStart', np.int32)
])

class GCMC:
    
    def __init__(self):

        self.ff_files = ['charmm36.ff/ffnonbonded.itp', 'charmm36.ff/nbfix.itp', 'charmm36.ff/silcs.itp']


        self.fragmentName = ['BENX', 'PRPX', 'DMEE', 'MEOH', 'FORM', 'IMIA', 'ACEY', 'MAMY', 'SOL']

        self.fragconc = [ 0.25,   0.25,   0.25,   0.25,   0.25,   0.25,   0.25,   0.25,  55.00]

        self.fragmuex = [-2.79, 1.46, -1.44, -5.36, -10.92, -14.18, -97.31, -68.49, -5.6]

        self.fragconf = [1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000]

        self.mctime = [1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 10.0]

        self.mcsteps = 10000

        self.seed = None

        self.attempt_prob_frag = [0.300, 0.300, 0.200, 0.200]
        
        self.attempt_prob_water = [0.250, 0.250, 0.250, 0.250]



        self.cutoff = 15.0

        self.fixCutoff = 6.0
        
        self.grid_dx = 1.0


        self.cavity_bias = True
        self.cavity_bias_factor = 1.0
        # self.cavity_bias_dx = 0.0

        # self.grid = np.zeros(1, dtype = np.int32)

        self.configurational_bias = True


        
        self.top_file = None
        self.pdb_file = None

        self.showInfo = False

        self.PME = False






        if os.path.exists('temp_link'):
            os.remove('temp_link')
        os.symlink(pkg_resources.resource_filename(__name__, 'charmm36.ff'), 'temp_link') # create a symbolic link to the force field directory
        os.rename('temp_link', 'charmm36.ff') # rename the symbolic link to force field directory


    def get_grid(self):
        print("Getting grid points...")

        x_values = [atom.x for atom in self.atoms]
        y_values = [atom.y for atom in self.atoms]
        z_values = [atom.z for atom in self.atoms]

        

        print("x direction: min =", min(x_values), "max =", max(x_values),end = '\t')
        print("y direction: min =", min(y_values), "max =", max(y_values),end = '\t')
        print("z direction: min =", min(z_values), "max =", max(z_values))

        x_center = (min(x_values) + max(x_values)) / 2.0
        y_center = (min(y_values) + max(y_values)) / 2.0
        z_center = (min(z_values) + max(z_values)) / 2.0

        print("x center =", x_center, end = '\t')
        print("y center =", y_center, end = '\t')
        print("z center =", z_center)

        self.startxyz = [x_center - self.cryst[0]/2.0, y_center - self.cryst[1]/2.0, z_center - self.cryst[2]/2.0]

        print("startxyz =", self.startxyz, end = '\t')

        self.endxyz = [x_center + self.cryst[0]/2.0, y_center + self.cryst[1]/2.0, z_center + self.cryst[2]/2.0]

        print("endxyz =", self.endxyz)

        # self.grid_n = [int((self.endxyz[0] - self.startxyz[0]) / self.grid_dx), int((self.endxyz[1] - self.startxyz[1]) / self.grid_dx), int((self.endxyz[2] - self.startxyz[2]) / self.grid_dx)]

        # print("grid_n =", self.grid_n)

        



        if self.cavity_bias:
            print("Using cavity bias: True")

            print("Grid dx =", self.grid_dx)
            self.grid_n = [int((self.endxyz[0] - self.startxyz[0]) / self.grid_dx), int((self.endxyz[1] - self.startxyz[1]) / self.grid_dx), int((self.endxyz[2] - self.startxyz[2]) / self.grid_dx)]
            print("grid_n =", self.grid_n)

            grid = {(x,y,z) for x in range(self.grid_n[0]) for y in range(self.grid_n[1]) for z in range(self.grid_n[2])}
            
            
            # change_dx = set()

            # for x in set(np.arange(0, self.cavity_bias_dx + self.grid_dx/2, self.grid_dx)) | set(np.arange(0,-self.cavity_bias_dx - self.grid_dx/2, -self.grid_dx)):
            #     for y in set(np.arange(0, self.cavity_bias_dx + self.grid_dx/2, self.grid_dx)) | set(np.arange(0,-self.cavity_bias_dx - self.grid_dx/2, -self.grid_dx)):
            #         for z in set(np.arange(0, self.cavity_bias_dx + self.grid_dx/2, self.grid_dx)) | set(np.arange(0,-self.cavity_bias_dx - self.grid_dx/2, -self.grid_dx)):
            #             change_dx.add((x, y, z))
            
            # print("change_dx =", change_dx)

            # for atom in self.atoms:
            #     for dx in change_dx:
            #         x = int((atom.x - self.startxyz[0] + dx[0]) / self.grid_dx) % self.grid_n[0]
            #         y = int((atom.y - self.startxyz[1] + dx[1]) / self.grid_dx) % self.grid_n[1]
            #         z = int((atom.z - self.startxyz[2] + dx[2]) / self.grid_dx) % self.grid_n[2]
            #         grid.discard(x + y * self.grid_n[0] + z * self.grid_n[0] * self.grid_n[1])
            for atom in self.atoms:
                x = int((atom.x - self.startxyz[0]) / self.grid_dx) % self.grid_n[0]
                y = int((atom.y - self.startxyz[1]) / self.grid_dx) % self.grid_n[1]
                z = int((atom.z - self.startxyz[2]) / self.grid_dx) % self.grid_n[2]
                grid.discard((x,y,z))
            
            

            # grid_sum = np.sum(self.grid)
            self.cavity_bias_factor = len(grid) / (self.grid_n[0] * self.grid_n[1] * self.grid_n[2])
            print("The number of grid points =", len(grid))
            print("cavity_bias_factor =", self.cavity_bias_factor)
            if self.cavity_bias_factor < 0.05:
                print("Error: cavity_bias_factor is too small, please decrease cavity_bias_dx")
                sys.exit()
        else:
            print("Using cavity bias: False")

            self.cavity_bias_factor = 1.0
            self.grid_dx = 1.0

            print("Grid dx =", self.grid_dx)
            self.grid_n = [int((self.endxyz[0] - self.startxyz[0]) / self.grid_dx), int((self.endxyz[1] - self.startxyz[1]) / self.grid_dx), int((self.endxyz[2] - self.startxyz[2]) / self.grid_dx)]
            print("grid_n =", self.grid_n)

            grid = {(x,y,z) for x in range(self.grid_n[0]) for y in range(self.grid_n[1]) for z in range(self.grid_n[2])}
            print("The number of grid points =", len(grid))
            print("cavity_bias_factor =", self.cavity_bias_factor)
            
        self.grid = np.zeros(len(grid) * 3, dtype = np.float32)
        for i, grid_point in enumerate(grid):
            x, y, z = grid_point
            self.grid[i*3] = x * self.grid_dx + self.startxyz[0]
            self.grid[i*3+1] = y * self.grid_dx + self.startxyz[1]
            self.grid[i*3+2] = z * self.grid_dx + self.startxyz[2]
        
        del grid
        
            # self.grid = np.zeros(self.grid_n[0] * self.grid_n[1] * self.grid_n[2], dtype = np.int32)
            # print("The number of grid points =", len(self.grid))
            # change_dx = [((i-1) * self.cavity_bias_dx,(j-1) * self.cavity_bias_dx, (k-1) * self.cavity_bias_dx) for i in range(3) for j in range(3) for k in range(3)]
            # change_dx = set(change_dx)
            # # print('change_dx =', change_dx)
            # for atom in self.atoms:
            #     for dx in change_dx:
            #         x = int((atom.x - self.startxyz[0] + dx[0]) / self.grid_dx) % self.grid_n[0]
            #         y = int((atom.y - self.startxyz[1] + dx[1]) / self.grid_dx) % self.grid_n[1]
            #         z = int((atom.z - self.startxyz[2] + dx[2]) / self.grid_dx) % self.grid_n[2]

            #         self.grid[x + y * self.grid_n[0] + z * self.grid_n[0] * self.grid_n[1]] += 1
            # grid_sum = sum(bool(x) for x in self.grid)
            # self.cavity_bias_factor = 1.0 - grid_sum / len(self.grid)
            # print("The number of grid points with atoms =", grid_sum)
            # print("The cavity bias factor =", self.cavity_bias_factor)
            # if self.cavity_bias_factor < 0.1 and self.cavity_bias_factor > 0.0:
            #     print("Error: The cavity bias factor is too small. Please decrease the cavity bias dx.")
            #     sys.exit(1) 


            
        
    
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

        self.volume = self.cryst[0] * self.cryst[1] * self.cryst[2]

        radius = self.cutoff / 2
        volume = (4 / 3) * np.pi * (radius ** 3)

        self.volume_n = int(self.volume / volume) + 1

        print(f"Total volume: {self.volume}")
        # print(f"Total volume: {self.volume}", end = '\t')
        # print(f"Single volume: {volume}", end = '\t')
        # print(f"volume_n: {self.volume_n}")

        
        








        
            
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
        
    def get_psf(self,psf_file):
        self.psf_file = psf_file

        print(f"Using psf file: {psf_file}")

        try:
            self.atom_top = protein_data.read_psf(psf_file)
        except:
            print(f"Error reading psf file: {psf_file}")
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

            print(f"Read solute: {frag} with {len(self.fragments[-1])} atoms")
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
        self.move_array = np.array(random.choices(range(len(self.fragmentName)), weights=self.mctime, k=self.mcsteps),dtype=np.int32)
        # self.move_array = np.random.choice(np.arange(len(self.fragmentName),dtype=np.int32), p=self.mctime, size=self.mcsteps)


        for i, n in enumerate(self.move_array):
            if self.fragmentName[n] == 'SOL':
                self.move_array[i] = self.move_array[i] * 4 + random.choices(range(len(self.attempt_prob_water)), weights=self.attempt_prob_water, k=1)[0]
            else:
                self.move_array[i] = self.move_array[i] * 4 + random.choices(range(len(self.attempt_prob_frag)), weights=self.attempt_prob_frag, k=1)[0]

        # for i, n in enumerate(self.move_array):
        #     if self.fragmentName[n] == 'SOL':
        #         self.move_array[i] = self.move_array[i] * 4 + np.random.choice(range(len(self.attempt_prob_water)), p=self.attempt_prob_water, size=1)[0]
        #     else:
        #         self.move_array[i] = self.move_array[i] * 4 + np.random.choice(range(len(self.attempt_prob_frag)), p=self.attempt_prob_frag, size=1)[0]


        self.move_array_n = [0 for i in range(len(self.attempt_prob_frag) * len(self.fragmentName))]
        self.move_array_frag = [0 for i in range(len(self.fragmentName))]

        for i in self.move_array:
            self.move_array_n[i] += 1
            self.move_array_frag[i//4] += 1
        
        movementName = ['Insert', 'Delete', 'Translate', 'Rotate']
        for i in range(len(self.fragmentName)):
            print(f"Solute {self.fragmentName[i]}\t Movement {self.move_array_frag[i]} times", end='\t')
            for j in range(len(self.attempt_prob_frag)):
                print(f"{movementName[j]} {self.move_array_n[i*4+j]} times", end='\t')
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

        for i in fragconf:
            if i <= 0:
                print("Error: fragconf number <= 0")
                sys.exit(1)

        if len(fragconf) == 1:
            fragconf = fragconf * len(self.fragmentName)
            self.fragconf = fragconf
        elif len(fragconf) == len(self.fragmentName):
            self.fragconf = fragconf
        else:
            print("Error: fragconf number not match")
            sys.exit(1)
        
        if all(x == 1 for x in self.fragconf):
            self.configurational_bias = False

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
        print("Solute Name: \t\t",'\t\t'.join(self.fragmentName))
        print("Solute Muex: \t\t",'\t\t'.join([str(i) for i in self.fragmuex]))
        print("Solute Conc: \t\t",'\t\t'.join([str(i) for i in self.fragconc]))
        print("Solute ConfB: \t",'\t\t'.join([str(i) for i in self.fragconf]))
        print("Solute mcTime: \t",'\t\t'.join([str(i) for i in self.mctime]))


        


    def get_simulation(self):

        
        self.get_fragment()
        self.get_forcefield()

        self.get_move()
        self.get_grid()

        self.show_parameters()

        
        
        if self.pdb_file is None:
            print("Error: pdb file not set")
            sys.exit(1)

        if self.top_file is None and self.psf_file is None:
            print("Error: top file or psf file not set")
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
        
        print(f"Solute atom type number: {len(self.atomtypes1)}")
        print(f"Total atom type number: {len(self.atomtypes2)}")

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
        
        # for i, frag in enumerate(self.fragments):
        #     print(f"Fragment {self.fragmentName[i]} number: {len(self.fraglist[i])}")


        self.set_fixCut()

        self.update_data()

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

    def update_data(self):


        
        
        self.fragmentInfo = np.empty(len(self.fragmentName), dtype=AtomArray_dtype)
        for i, frag in enumerate(self.fragments):
            # self.fragmentInfo[i]['name'] = self.fragmentName[i]
            self.fragmentInfo[i]['name'] = np.array(self.fragmentName[i].ljust(4)[:4], dtype='S4')

            
            self.fragmentInfo[i]['muex'] = self.fragmuex[i]
            self.fragmentInfo[i]['conc'] = self.fragconc[i]
            self.fragmentInfo[i]['confBias'] = self.fragconf[i]
            self.fragmentInfo[i]['mcTime'] = self.mctime[i]

            self.fragmentInfo[i]['totalNum'] = len(self.fraglist[i])

            maxNum1 = self.fragconc[i] * self.volume * MOLES_TO_MOLECULES + self.move_array_n[i*4]
            maxNum2 = len(self.fraglist[i]) + self.move_array_n[i*4] 

            self.fragmentInfo[i]['maxNum'] = max(maxNum1, maxNum2)

            print(f"Solute {self.fragmentName[i]}: Total number: {self.fragmentInfo[i]['totalNum']}, Max number: {self.fragmentInfo[i]['maxNum']}")

            self.fragmentInfo[i]['num_atoms'] = len(frag)

            atom_center = np.zeros(3)
            for atom in frag:
                atom_center += np.array([atom.x, atom.y, atom.z])
            atom_center /= len(frag)
            # print(f"Fragment {self.fragmentName[i]}: Center of conf: {atom_center}", end='\t')
            
            for atom in frag:
                atom.x -= atom_center[0]
                atom.y -= atom_center[1]
                atom.z -= atom_center[2]
            # atom_center = np.zeros(3)

            for j, atom in enumerate(frag):
                self.fragmentInfo[i]['atoms'][j]['position'][0] = atom.x
                self.fragmentInfo[i]['atoms'][j]['position'][1] = atom.y
                self.fragmentInfo[i]['atoms'][j]['position'][2] = atom.z
                self.fragmentInfo[i]['atoms'][j]['charge'] = atom.charge
                self.fragmentInfo[i]['atoms'][j]['type'] = atom.typeNum

                # atom_center += np.array([atom.x, atom.y, atom.z])
        if len(self.fix_atoms) == 0:
            TotalResidueNum = sum([self.fragmentInfo[i]['maxNum'] for i in range(len(self.fragmentName))])
        else:
            TotalResidueNum = self.fix_atoms[-1].sequence2 + 1 + sum([self.fragmentInfo[i]['maxNum'] for i in range(len(self.fragmentName))])
        TotalAtomNum = len(self.fix_atoms) + sum([self.fragmentInfo[i]['maxNum'] * self.fragmentInfo[i]['num_atoms'] for i in range(len(self.fragmentName))]) 

        self.residueInfo = np.empty(TotalResidueNum, dtype=Residue_dtype)
        self.atomInfo = np.empty(TotalAtomNum, dtype=Atom_dtype)

        n = -1
        for i, atom in enumerate(self.fix_atoms):

            sequence = atom.sequence2 
            if sequence != n:
                n = sequence
                self.residueInfo[sequence]['atomNum'] = 0
                self.residueInfo[sequence]['atomStart'] = i
                self.residueInfo[sequence]['position'][0] = 0
                self.residueInfo[sequence]['position'][1] = 0
                self.residueInfo[sequence]['position'][2] = 0
            
            self.residueInfo[sequence]['atomNum'] += 1
            self.residueInfo[sequence]['position'][0] += atom.x
            self.residueInfo[sequence]['position'][1] += atom.y
            self.residueInfo[sequence]['position'][2] += atom.z

            self.atomInfo[i]['position'][0] = atom.x
            self.atomInfo[i]['position'][1] = atom.y
            self.atomInfo[i]['position'][2] = atom.z
            self.atomInfo[i]['charge'] = atom.charge
            self.atomInfo[i]['type'] = atom.typeNum
        if len(self.fix_atoms) == 0:
            sequence = -1
        for i in range(sequence + 1):
            self.residueInfo[i]['position'][0] /= self.residueInfo[i]['atomNum']
            self.residueInfo[i]['position'][1] /= self.residueInfo[i]['atomNum']
            self.residueInfo[i]['position'][2] /= self.residueInfo[i]['atomNum']
        
        for i, frag in enumerate(self.fragments):
            self.fragmentInfo[i]['startRes'] = sequence + 1
            for j in range(self.fragmentInfo[i]['maxNum']):
                sequence += 1
                self.residueInfo[sequence]['atomNum'] = self.fragmentInfo[i]['num_atoms']
                self.residueInfo[sequence]['atomStart'] = self.residueInfo[sequence - 1]['atomStart'] + self.residueInfo[sequence - 1]['atomNum']
                self.residueInfo[sequence]['position'][0] = 0
                self.residueInfo[sequence]['position'][1] = 0
                self.residueInfo[sequence]['position'][2] = 0
                if j < self.fragmentInfo[i]['totalNum']:
                    for k, atom in enumerate(self.fraglist[i][j]):
                        self.residueInfo[sequence]['position'][0] += atom.x
                        self.residueInfo[sequence]['position'][1] += atom.y
                        self.residueInfo[sequence]['position'][2] += atom.z

                        self.atomInfo[self.residueInfo[sequence]['atomStart'] + k]['position'][0] = atom.x
                        self.atomInfo[self.residueInfo[sequence]['atomStart'] + k]['position'][1] = atom.y
                        self.atomInfo[self.residueInfo[sequence]['atomStart'] + k]['position'][2] = atom.z
                        self.atomInfo[self.residueInfo[sequence]['atomStart'] + k]['charge'] = atom.charge
                        self.atomInfo[self.residueInfo[sequence]['atomStart'] + k]['type'] = atom.typeNum
                    
                    self.residueInfo[sequence]['position'][0] /= self.residueInfo[sequence]['atomNum']
                    self.residueInfo[sequence]['position'][1] /= self.residueInfo[sequence]['atomNum']
                    self.residueInfo[sequence]['position'][2] /= self.residueInfo[sequence]['atomNum']
                

                

                # self.residueInfo[sequence]['position'][0] /= self.residueInfo[sequence]['atomNum']
                # self.residueInfo[sequence]['position'][1] /= self.residueInfo[sequence]['atomNum']
                # self.residueInfo[sequence]['position'][2] /= self.residueInfo[sequence]['atomNum']

        
        self.ff = np.empty(len(self.ff_pairs) * 2, dtype=np.float32)
        for i, pair in enumerate(self.ff_pairs):
            self.ff[i * 2] = pair[0]
            self.ff[i * 2 + 1] = pair[1]
        
        
        self.SimInfo = np.empty(1, dtype=Info_dtype)


        self.SimInfo[0]['ffXNum'] = len(self.atomtypes1)
        self.SimInfo[0]['ffYNum'] = len(self.atomtypes2)

        self.SimInfo[0]['mcsteps'] = self.mcsteps
        self.SimInfo[0]['cutoff'] = self.cutoff
        self.SimInfo[0]['grid_dx'] = self.grid_dx
        self.SimInfo[0]['startxyz'][0] = self.startxyz[0]
        self.SimInfo[0]['startxyz'][1] = self.startxyz[1]
        self.SimInfo[0]['startxyz'][2] = self.startxyz[2]
        self.SimInfo[0]['cryst'][0] = self.cryst[0]
        self.SimInfo[0]['cryst'][1] = self.cryst[1]
        self.SimInfo[0]['cryst'][2] = self.cryst[2]

        self.SimInfo[0]['cavityFactor'] = self.cavity_bias_factor

        self.SimInfo[0]['fragTypeNum'] = len(self.fragmentName)

        self.SimInfo[0]['totalGridNum'] = len(self.grid) // 3
        self.SimInfo[0]['totalResNum'] = TotalResidueNum
        self.SimInfo[0]['totalAtomNum'] = TotalAtomNum

        self.SimInfo[0]['showInfo'] = self.showInfo

        self.SimInfo[0]['PME'] = self.PME

        self.SimInfo[0]['seed'] = random.randint(0, (2**32)-1)

        # 
        # print('showInfo', self.SimInfo[0]['showInfo'])

        # n = 0
        # for frag in self.fragmentInfo:
        #     for i in range(frag['num_atoms']):
        #         for atom in self.fix_atoms:
        #             n +=1
        #             if n % 100 == 0:
        #                 type1 = frag['atoms'][i]['type']
        #                 type2 = atom.typeNum
        #                 # print(type1,type2,self.atomtypes1[type1])
        #                 print(type1, type2, self.atomtypes1[type1], self.atomtypes2[type2], self.ff[(type1 * self.SimInfo[0]['ffXNum'] + type2)* 2], self.ff[(type1 * self.SimInfo[0]['ffXNum'] + type2)* 2 + 1])
                

        # for i in self.atomtypes1:
        #     print(i, end=' ')
        # print()
        # for i in self.atomtypes2:
        #     print(i, end=' ')
        # print()

        


            # atom_center /= len(frag)
            # print(f"Fragment {self.fragmentName[i]}: The final center of conf: {atom_center}")
        
    def run(self):
        self.starttime = time.time()
        print('Start GPU simulation...')

        runGCMC(self.SimInfo, self.fragmentInfo, self.residueInfo, self.atomInfo, self.grid, self.ff, self.move_array)

        self.endtime = time.time()
        print('End GPU simulation...')
        print('GPU simulation time: %s s' % (self.endtime - self.starttime))

def main():
    # file_output = open('Analyze_output.txt', 'w')
    # original_output = sys.stdout
    # sys.stdout = Tee(sys.stdout, file_output)

    starttime = time.time()
    print('Start GCMC simulation at %s...' % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    

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
        "-s",
        "--psf-file",
        dest="psf_file",
        required=False,
        help="The file .psf for GCMC",
        metavar="file.psf",
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
        help="The value of solute muex(splice by , with no space), if the first value is negative, then follow the -u or --fragmuex without space",
        metavar="muex1,muex2,...",
        type=str,
    ) 
    parser.add_argument(
        "-f",
        "--fragconf",
        dest="fragconf",
        required=False,
        help="The value of solute conf(splice by , with no space). Or only one value for all solutes",
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
        help="The mctime of solutes(splice by , with no space)",
        metavar="mctime1,mctime2,...",
        type=str,
    )
    parser.add_argument(
        "-c",
        "--fragconc",
        dest="fragconc",
        required=False,
        help="The value of solute concentration(splice by , with no space)",
        metavar="conc1,conc2,...",
        type=str,
    )

    parser.add_argument(
        "-y",
        "--cavitybias-dx",
        dest="cavity_bias_dx",
        required=False,
        help="The value of cavity bias dx(if dx <= 0, then no cavity bias)",
        metavar="cavity_bias_dx",
        type=float,
    )

    parser.add_argument(
        "-e",
        "--seed",
        dest="seed",
        required=False,
        help="The seed of random number",
        metavar="seed",
        type=int,
    )

    parser.add_argument(
        "-P",
        "--PME",
        dest="PME",
        required=False,
        help="Enable PME(Default: Disable)",
        action='store_true',
    )


    parser.add_argument(
        "-w",
        "--show-info",
        dest="show_info",
        required=False,
        help="Show the information of solutes",
        action='store_true',
    )


    args = parser.parse_args()
    
    gcmc = GCMC()

    pdb_file = args.pdb_file
    top_file = args.top_file
    out_file = args.out_file
    psf_file = args.psf_file

    if out_file is not None:
        file_output = open(out_file, 'w')
        original_output = sys.stdout
        sys.stdout = Tee(sys.stdout, file_output)
        print(f"Using output file: {out_file}")

    if args.fragmuex is not None:
        fragmuex = args.fragmuex
        fragmuex = fragmuex.split(',')
        gcmc.get_fragmuex(fragmuex)
        print(f"Using solute muex: {fragmuex}")    

    if args.fragconf is not None:
        fragconf = args.fragconf
        fragconf = fragconf.split(',')
        gcmc.get_fragconf(fragconf)
        print(f"Using solute conf: {fragconf}")

    if args.mcsteps is not None:
        gcmc.mcsteps = args.mcsteps
        print(f"Using MC steps: {args.mcsteps}")

    if args.mctime is not None:
        mctime = args.mctime
        mctime = mctime.split(',')
        gcmc.get_mctime(mctime)
        print(f"Using solute mctime: {mctime}")

    if args.fragconc is not None:
        fragconc = args.fragconc
        fragconc = fragconc.split(',')
        gcmc.get_fragconc(fragconc)
        print(f"Using solute concentration: {fragconc}")
    
    if args.cavity_bias_dx is not None:
        gcmc.grid_dx = args.cavity_bias_dx
        if gcmc.grid_dx <= 0:
            gcmc.cavity_bias = False
            print(f"No cavity bias")
        else:
            print(f"Using cavity bias dx: {args.cavity_bias_dx}")

    if args.seed is not None:
        gcmc.seed = args.seed
        random.seed(gcmc.seed)
        print(f"Using seed: {args.seed}")

    if args.PME:
        gcmc.PME = True
        print(f"Using PME")
    
    if args.show_info:
        gcmc.showInfo = True
        print(f"Showing the information of solutes")

    gcmc.get_pdb(pdb_file)

        
    
    



    if top_file is not None:
        gcmc.get_top(top_file)
    
    if psf_file is not None:
        gcmc.get_psf(psf_file)


    
    gcmc.get_simulation()


    endtime = time.time()
    print(f"Python time used: {endtime - starttime} s")

    gcmc.run()

    
    print('GCMC simulation finished at %s...' % time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

    endtime = time.time()
    print(f"Time used: {endtime - starttime} s")


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