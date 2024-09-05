"""
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu
"""

import sys
import numpy as np
import random
from .values import *
import copy
import re


class GCMCDataset:
    
    def get_grid(self):
        # Generate grid points for the simulation
        print("Getting grid points...")

        x_values = [atom.x for atom in self.atoms]
        y_values = [atom.y for atom in self.atoms]
        z_values = [atom.z for atom in self.atoms]

        # Print grid dimensions
        print("x direction: min =", min(x_values), "max =", max(x_values),end = '\t')
        print("y direction: min =", min(y_values), "max =", max(y_values),end = '\t')
        print("z direction: min =", min(z_values), "max =", max(z_values))

        # Calculate center of the grid
        x_center = (min(x_values) + max(x_values)) / 2.0
        y_center = (min(y_values) + max(y_values)) / 2.0
        z_center = (min(z_values) + max(z_values)) / 2.0

        print("x center =", x_center, end = '\t')
        print("y center =", y_center, end = '\t')
        print("z center =", z_center)

        # Set start and end coordinates of the grid
        self.startxyz = [x_center - self.cryst[0]/2.0, y_center - self.cryst[1]/2.0, z_center - self.cryst[2]/2.0]
        print("startxyz =", self.startxyz, end = '\t')
        self.endxyz = [x_center + self.cryst[0]/2.0, y_center + self.cryst[1]/2.0, z_center + self.cryst[2]/2.0]
        print("endxyz =", self.endxyz)

        if self.cavity_bias:
            # Generate grid with cavity bias
            print("Using cavity bias: True")
            print("Grid dx =", self.grid_dx)
            self.grid_n = [int((self.endxyz[0] - self.startxyz[0]) / self.grid_dx), int((self.endxyz[1] - self.startxyz[1]) / self.grid_dx), int((self.endxyz[2] - self.startxyz[2]) / self.grid_dx)]
            print("grid_n =", self.grid_n)

            grid = {(x,y,z) for x in range(self.grid_n[0]) for y in range(self.grid_n[1]) for z in range(self.grid_n[2])}
            
            # Remove grid points occupied by atoms
            for atom in self.atoms:
                x = int((atom.x - self.startxyz[0]) / self.grid_dx) % self.grid_n[0]
                y = int((atom.y - self.startxyz[1]) / self.grid_dx) % self.grid_n[1]
                z = int((atom.z - self.startxyz[2]) / self.grid_dx) % self.grid_n[2]
                grid.discard((x,y,z))
            
            # Calculate cavity bias factor
            self.cavity_bias_factor = len(grid) / (self.grid_n[0] * self.grid_n[1] * self.grid_n[2])
            print("The number of grid points =", len(grid))
            print("cavity_bias_factor =", self.cavity_bias_factor)
            if self.cavity_bias_factor < 0.05:
                print("Error: cavity_bias_factor is too small, please decrease cavity_bias_dx")
                sys.exit()
        else:
            # Generate uniform grid without cavity bias
            print("Using cavity bias: False")
            self.cavity_bias_factor = 1.0
            self.grid_dx = 1.0
            print("Grid dx =", self.grid_dx)
            self.grid_n = [int((self.endxyz[0] - self.startxyz[0]) / self.grid_dx), int((self.endxyz[1] - self.startxyz[1]) / self.grid_dx), int((self.endxyz[2] - self.startxyz[2]) / self.grid_dx)]
            print("grid_n =", self.grid_n)
            grid = {(x,y,z) for x in range(self.grid_n[0]) for y in range(self.grid_n[1]) for z in range(self.grid_n[2])}
            print("The number of grid points =", len(grid))
            print("cavity_bias_factor =", self.cavity_bias_factor)
            
        # Create grid array
        self.grid = np.zeros(len(grid) * 3, dtype = np.float32)
        for i, grid_point in enumerate(grid):
            x, y, z = grid_point
            self.grid[i*3] = x * self.grid_dx + self.startxyz[0]
            self.grid[i*3+1] = y * self.grid_dx + self.startxyz[1]
            self.grid[i*3+2] = z * self.grid_dx + self.startxyz[2]
        
        del grid
        
    def get_move(self):
        # Generate move array for Monte Carlo simulation
        self.moveArray = np.array(random.choices(range(len(self.fragmentName)), weights=self.mctime, k=self.mcsteps),dtype=np.int32)

        # Assign move types based on fragment type
        for i, n in enumerate(self.moveArray):
            if self.fragmentName[n] == 'SOL':
                self.moveArray[i] = self.moveArray[i] * 4 + random.choices(range(len(self.attempt_prob_water)), weights=self.attempt_prob_water, k=1)[0]
            else:
                self.moveArray[i] = self.moveArray[i] * 4 + random.choices(range(len(self.attempt_prob_frag)), weights=self.attempt_prob_frag, k=1)[0]

        # Calculate move statistics
        self.moveArray_n = [0 for i in range(len(self.attempt_prob_frag) * len(self.fragmentName))]
        self.moveArray_frag = [0 for i in range(len(self.fragmentName))]

        for i in self.moveArray:
            self.moveArray_n[i] += 1
            self.moveArray_frag[i//4] += 1
        
        # Print move statistics
        movementName = ['Insert', 'Delete', 'Translate', 'Rotate']
        for i in range(len(self.fragmentName)):
            print(f"Solute {self.fragmentName[i]}\t Movement {self.moveArray_frag[i]} times", end='\t')
            for j in range(len(self.attempt_prob_frag)):
                print(f"{movementName[j]} {self.moveArray_n[i*4+j]} times", end='\t')
            print()

    def set_fixCut(self):
        # Set fixed atoms cutoff
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
        # Update simulation data structures
        
        # Update fragment information
        self.fragmentInfo = np.empty(len(self.fragmentName), dtype=AtomArray_dtype)
        for i, frag in enumerate(self.fragments):
            self.fragmentInfo[i]['name'] = np.array(self.fragmentName[i].ljust(4)[:4], dtype='S4')
            self.fragmentInfo[i]['muex'] = self.fragmuex[i]
            self.fragmentInfo[i]['conc'] = self.fragconc[i]
            self.fragmentInfo[i]['confBias'] = self.fragconf[i]
            self.fragmentInfo[i]['mcTime'] = self.mctime[i]
            self.fragmentInfo[i]['totalNum'] = len(self.fraglist[i])

            # Calculate maximum number of fragments
            maxNum1 = self.fragconc[i] * self.volume * MOLES_TO_MOLECULES * 2 + self.moveArray_n[i*4]
            maxNum2 = len(self.fraglist[i]) + self.moveArray_n[i*4] 
            self.fragmentInfo[i]['maxNum'] = max(maxNum1, maxNum2)

            print(f"Solute {self.fragmentName[i]}: Total number: {self.fragmentInfo[i]['totalNum']}, Max number: {self.fragmentInfo[i]['maxNum']}")

            self.fragmentInfo[i]['num_atoms'] = len(frag)

            # Center fragment atoms
            atom_center = np.zeros(3)
            for atom in frag:
                atom_center += np.array([atom.x, atom.y, atom.z])
            atom_center /= len(frag)
            
            for atom in frag:
                atom.x -= atom_center[0]
                atom.y -= atom_center[1]
                atom.z -= atom_center[2]

            for j, atom in enumerate(frag):
                self.fragmentInfo[i]['atoms'][j]['position'][0] = atom.x
                self.fragmentInfo[i]['atoms'][j]['position'][1] = atom.y
                self.fragmentInfo[i]['atoms'][j]['position'][2] = atom.z
                self.fragmentInfo[i]['atoms'][j]['charge'] = atom.charge
                self.fragmentInfo[i]['atoms'][j]['type'] = atom.typeNum

        # Calculate total residue and atom numbers
        if len(self.fix_atoms) == 0:
            TotalResidueNum = sum([self.fragmentInfo[i]['maxNum'] for i in range(len(self.fragmentName))])
        else:
            TotalResidueNum = self.fix_atoms[-1].sequence2 + 1 + sum([self.fragmentInfo[i]['maxNum'] for i in range(len(self.fragmentName))])
        TotalAtomNum = len(self.fix_atoms) + sum([self.fragmentInfo[i]['maxNum'] * self.fragmentInfo[i]['num_atoms'] for i in range(len(self.fragmentName))]) 

        # Initialize residue and atom information arrays
        self.residueInfo = np.empty(TotalResidueNum, dtype=Residue_dtype)
        self.atomInfo = np.empty(TotalAtomNum, dtype=Atom_dtype)

        # Update fixed atoms information
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
        
        # Update fragment atoms information
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
        
        # Update force field information
        self.ff = np.empty(len(self.ff_pairs) * 2, dtype=np.float32)
        for i, pair in enumerate(self.ff_pairs):
            self.ff[i * 2] = pair[0]
            self.ff[i * 2 + 1] = pair[1]
        
        # Update simulation information
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

    def get_data(self):
        # Prepare output data for the simulation

        print("Getting data...")
        
        for i in range(len(self.fragmentName)):
            print(f"Solute %s: Total number: %d" % (self.fragmentName[i], self.fragmentInfo[i]['totalNum']))

        # Prepare PDB string
        s = 'CRYST1  %.3f  %.3f  %.3f  90.00  90.00  90.00 P 1           1\n' % (self.cryst[0], self.cryst[1], self.cryst[2])

        for atom in self.fix_atoms:
            s += atom.s
            s += '\n'

        try:
            atomNum = int(self.fix_atoms[-1].serial) + 1
        except:
            atomNum = 1

        try:
            residueNum = int(self.fix_atoms[-1].sequence) + 1
        except:
            residueNum = 1

        for i in range(len(self.fragmentName)):
            for j in range(self.fragmentInfo[i]['totalNum']):
                resNum = self.fragmentInfo[i]['startRes'] + j
                res = self.residueInfo[resNum]
                for k in range(res['atomNum']):
                    atom = self.atomInfo[res['atomStart'] + k]
                    s += 'ATOM  %5d  %-4s%-4s %4d    %8.3f%8.3f%8.3f  1.00  0.00\n' % ((atomNum - 1) % 99999 + 1, self.fragments[i][k].name[:4], self.fragmentName[i][:4], (residueNum - 1) % 9999 + 1, atom['position'][0], atom['position'][1], atom['position'][2])
                    atomNum += 1
                residueNum += 1
        s += 'END\n'

        self.PDBString = s

        # Update topology file if needed
        if self.topologyType == 'top':
            s = copy.deepcopy(self.TOPString)
            try:
                pattern = r'\[\s*molecules\s*\]'
                parts = re.split(pattern, s)
                parts[1] = parts[1].strip()+'\n'
                p = re.findall(r'(\S+?)\s+(\d+)', parts[1])
                d = {i[0]:i[1] for i in p}

                for i in range(len(self.fragmentName)):
                    if self.fragmentName[i] in d:
                        ni = self.fragmentName[i]
                        parts[1] = re.sub('%s\\s+%s' % (ni,d[ni]), '%-s\\t%d' % (self.fragmentName[i], self.fragmentInfo[i]['totalNum']), parts[1])  
                    else:
                        parts[1] += '%-s\t%d\n' % (self.fragmentName[i], self.fragmentInfo[i]['totalNum'])
                
                self.TOPString = parts[0] + '\n[ molecules ]\n' + parts[1]

            except:
                print("Error: writing top file error")
                sys.exit(1)