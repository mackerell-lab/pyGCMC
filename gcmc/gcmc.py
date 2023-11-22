"""

    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Abhishek A. Kognole, 
    	Aoxiang Tao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu

"""

# from .main import main
# from .mainOld import main as mainOld

# from .packages import *
# from .values import *
from .base import GCMCBase
from .files import GCMCFiles
from .parameters import GCMCParameters
from .simulation import GCMCSimulation
from .dataset import GCMCDataset

import time
from .gpu import runGCMC


class GCMC(GCMCBase,GCMCFiles,GCMCParameters,GCMCDataset,GCMCSimulation):
    
    def __init__(self):

        GCMCBase.__init__(self)

    def run(self):
        self.starttime = time.time()
        print('Start GPU simulation...')

        runGCMC(self.SimInfo, self.fragmentInfo, self.residueInfo, self.atomInfo, self.grid, self.ff, self.move_array)

        self.endtime = time.time()
        print('End GPU simulation...')
        print('GPU simulation time: %s s' % (self.endtime - self.starttime))




