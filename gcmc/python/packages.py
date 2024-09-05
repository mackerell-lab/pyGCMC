"""
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu
"""


import os
import sys
import time
import random
import argparse
import numpy as np
import pkg_resources

from . import protein_data
from ..gpu import runGCMC

# https://stackoverflow.com/questions/616645/how-do-i-duplicate-sys-stdout-to-a-log-file-in-python
# Class to duplicate stdout to a log file
class Tee(object):
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()  # For immediate output visibility

    def flush(self):
        for f in self.files:
            f.flush()
