"""
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu
"""

import sys
import numpy as np
from . import protein_data

class GCMCFiles:

    # Load and parse configuration parameters
    def load_parameters(self, filename):
        self.config_dict = {}
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if line:
                    line = line.split('#', 1)[0].split('//', 1)[0].split('!', 1)[0].strip()
                    if line:
                        key, value = line.split(':', 1)
                        key, value = key.strip(), value.strip()
                        self.config_dict.setdefault(key, []).append(value)
        
        for i, v in self.config_dict.items():
            print(i, v)

    # Convert string to boolean
    def bool_parameters(self, value):
        return value.lower() in ('yes', 'true', 'on', '1')

    # Read and process PDB file
    def get_pdb(self, pdb_file):
        self.pdb_file = pdb_file
        print(f"Using pdb file: {pdb_file}")

        try:
            self.cryst, self.atoms, self.PDBString = protein_data.read_pdb(pdb_file)
        except Exception as e:
            print(f"Error reading pdb file: {pdb_file}")
            print(f"Error message: {str(e)}")
            sys.exit(1)

        print(f"pdb atom number: {len(self.atoms)}")   
        print(f"pdb cryst: {self.cryst}")

        self.volume = self.cryst[0] * self.cryst[1] * self.cryst[2]
        radius = self.cutoff / 2
        volume = (4 / 3) * np.pi * (radius ** 3)
        self.volume_n = int(self.volume / volume) + 1

        print(f"Total volume: {self.volume}")

    # Read and process TOP file
    def get_top(self, top_file):
        self.topologyType = 'top'
        self.top_file = top_file
        print(f"Using top file: {top_file}")

        try:
            self.atom_top, self.TOPString = protein_data.read_top(top_file)
        except:
            print(f"Error reading top file: {top_file}")
            sys.exit(1)

        print(f"top atom number: {len(self.atom_top)}")

        if len(self.atoms) != len(self.atom_top):
            print(f"Error: pdb atom number {len(self.atoms)} != top atom number {len(self.atom_top)}")
            sys.exit(1)
    
        self._update_atom_info()

    # Read and process PSF file
    def get_psf(self, psf_file):
        self.topologyType = 'psf'
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
    
        self._update_atom_info()

    # Update atom information based on topology
    def _update_atom_info(self):
        for i, atom in enumerate(self.atoms):
            if atom.residue[:3] != self.atom_top[i][2][:3]:
                print(f"Error: pdb atom {i+1} residue {atom.residue} != top atom residue {self.atom_top[i][2]}")
                sys.exit(1)
            atom.residue = self.atom_top[i][2]
            atom.type = self.atom_top[i][0]
            atom.charge = self.atom_top[i][4]
            atom.nameTop = self.atom_top[i][3]

    # Read and process fragment files
    def get_fragment(self):
        self.fragments = []

        for frag in self.fragmentName:
            try:
                self.fragments += [protein_data.read_pdb(f"charmm36.ff/mol/{frag.lower()}.pdb")[1]]
                atom_top = protein_data.read_itp(f"charmm36.ff/mol/{frag.lower()}.itp")
                for i, atom in enumerate(self.fragments[-1]):
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

    # Read and process force field files
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

