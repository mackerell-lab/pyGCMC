"""
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu
"""

import numpy as np

# Conversion constant from moles to molecules
MOLES_TO_MOLECULES = 0.0006023

try:
    from .gcmcH import *
except ImportError:
    print("Error: Unable to import gcmcH.py. Please ensure the file exists and is in the correct path.")

