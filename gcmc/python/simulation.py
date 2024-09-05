"""
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu
"""

import sys
import numpy as np

class GCMCSimulation:
    
    def get_simulation(self):
        # Initialize simulation components
        self.get_fragment()
        self.get_forcefield()
        self.get_move()
        self.get_grid()
        self.show_parameters()

        # Check for required files
        self.check_required_files()
        
        # Process atom types
        self.process_atom_types()
        
        # Generate force field pairs
        self.generate_ff_pairs()
        
        # Assign type numbers to atoms
        self.assign_type_numbers()
        
        # Categorize atoms as fixed or relaxed
        self.categorize_atoms()
        
        # Set fixed atom cutoff
        self.set_fixCut()
        
        # Update simulation data
        self.update_data()

    def check_required_files(self):
        if self.pdb_file is None:
            print("Error: pdb file not set")
            sys.exit(1)
        if self.top_file is None and self.psf_file is None:
            print("Error: top file or psf file not set")
            sys.exit(1)

    def process_atom_types(self):
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

    def generate_ff_pairs(self):
        self.ff_pairs = []
        for type1 in self.atomtypes1:
            for type2 in self.atomtypes2:
                ff_pair = self.get_ff_pair(type1, type2)
                self.ff_pairs.append(ff_pair)
        
        print(f"Total FF pair number: {len(self.ff_pairs)}")

    def get_ff_pair(self, type1, type2):
        if (type1, type2) in self.nbfix_dict:
            return self.nbfix_dict[(type1, type2)]
        elif (type2, type1) in self.nbfix_dict:
            return self.nbfix_dict[(type2, type1)]
        else:
            sigma1, epsilon1 = self.nb_dict[type1]
            sigma2, epsilon2 = self.nb_dict[type2]
            sigma = (sigma1 + sigma2) / 2
            epsilon = np.sqrt(epsilon1 * epsilon2)
            return [sigma, epsilon]

    def assign_type_numbers(self):
        for atom in self.atoms:
            atom.typeNum = self.atomtypes2.index(atom.type)
        
        for frag in self.fragments:
            for atom in frag:
                atom.typeNum = self.atomtypes2.index(atom.type)

    def categorize_atoms(self):
        self.fix_atoms = []
        relax_atoms = [[] for _ in range(len(self.fragments))]
        for atom in self.atoms:
            if atom.residue in self.fragmentName:
                relax_atoms[self.fragmentName.index(atom.residue)].append(atom)
            else:
                self.fix_atoms.append(atom)
        
        self.fraglist = [
            [relax_atoms[i][j:j + len(self.fragments[i])] 
             for j in range(0, len(relax_atoms[i]), len(self.fragments[i]))]
            for i in range(len(self.fragments))
        ]
        
        print(f"Total fixed atom number: {len(self.fix_atoms)}")
        print(f"Total relax atom number: {len(self.atoms) - len(self.fix_atoms)}")
