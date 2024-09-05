"""
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu
"""


import os
import importlib.util

# Check if gcmcH.py exists
gcmch_path = os.path.join(os.path.dirname(__file__), 'python', 'gcmcH.py')
if not os.path.exists(gcmch_path):
    # If it doesn't exist, run update_values.py to generate it
    update_values_path = os.path.join(os.path.dirname(__file__), 'scripts', 'update_values.py')
    spec = importlib.util.spec_from_file_location("update_values", update_values_path)
    update_values = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(update_values)
    
    # Run the update_values_py function
    header_file = os.path.join(os.path.dirname(__file__), 'cpp', 'gcmc.h')
    update_values.update_values_py(header_file, gcmch_path)
    print(f"Generated {gcmch_path}")

from .python import GCMC
from .python import packages

