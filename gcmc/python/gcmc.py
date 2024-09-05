"""
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu
"""

import time
from .base import GCMCBase
from .files import GCMCFiles
from .parameters import GCMCParameters
from .simulation import GCMCSimulation
from .dataset import GCMCDataset
from ..gpu import runGCMC

# Main GCMC class integrating all modules
class GCMC(GCMCBase, GCMCFiles, GCMCParameters, GCMCDataset, GCMCSimulation):
    """
    GCMC class that integrates base, files, parameters, dataset, and simulation modules
    to perform a Grand Canonical Monte Carlo simulation.
    """

    # Initialize GCMC object
    def __init__(self):
        super().__init__()

    # Run GCMC simulation
    def run(self):
        self.startTime = time.time()
        print('Start GPU simulation...')

        runGCMC(
            self.SimInfo,
            self.fragmentInfo,
            self.residueInfo,
            self.atomInfo,
            self.grid,
            self.ff,
            self.moveArray
        )

        self.endTime = time.time()
        print('End GPU simulation...')
        print(f'GPU simulation time: {self.endTime - self.startTime:.2f} s')
