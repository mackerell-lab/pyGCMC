"""
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu
"""


import os
import pkg_resources

class GCMCBase:
    def __init__(self):
        
        # Force field files
        self.ff_files = ['charmm36.ff/ffnonbonded.itp', 'charmm36.ff/nbfix.itp', 'charmm36.ff/silcs.itp']
        
        # Fragment properties
        self.fragmentName = ['BENX', 'PRPX', 'DMEE', 'MEOH', 'FORM', 'IMIA', 'ACEY', 'MAMY', 'SOL']
        self.fragconc = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 55.00]
        self.fragmuex = [-2.79, 1.46, -1.44, -5.36, -10.92, -14.18, -97.31, -68.49, -5.6]
        self.fragconf = [1000] * 9
        self.mctime = [1.000] * 8 + [10.0]

        # Monte Carlo parameters
        self.mcsteps = 10000
        self.seed = None
        self.attempt_prob_frag = [0.300, 0.300, 0.200, 0.200]
        self.attempt_prob_water = [0.250, 0.250, 0.250, 0.250]

        # Simulation parameters
        self.cutoff = 15.0
        self.fixCutoff = 6.0
        self.grid_dx = 1.0
        self.cavity_bias = True
        self.cavity_bias_factor = 1.0
        self.configurational_bias = True

        # Input files
        self.top_file = None
        self.pdb_file = None

        # Simulation options
        self.showInfo = False
        self.PME = False

        # Set up force field directory
        self._setup_ff_directory()

    # Set up force field directory
    def _setup_ff_directory(self):
        if os.path.exists('temp_link'):
            os.remove('temp_link')
        os.symlink(pkg_resources.resource_filename(__name__, 'resources/charmm36.ff'), 'temp_link')
        os.rename('temp_link', 'charmm36.ff')


