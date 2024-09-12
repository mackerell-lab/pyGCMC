"""
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu
"""


import os
import importlib.util
import zipfile
# import pkg_resources




# Modify the way resource file paths are located
toppar_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources', 'toppar')
resources_zip = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources', 'resources.zip')
resources_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'resources')

# Check if CHARMM force field parameters exist and extract if necessary
if not os.path.exists(toppar_dir) or not os.path.exists(os.path.join(toppar_dir, 'silcs.prm')) or os.path.getsize(os.path.join(toppar_dir, 'silcs.prm')) < 10:
    print("CHARMM force field parameters not found.")
    print("Extracting CHARMM force field parameters...")
    if os.path.exists(resources_zip):
        with zipfile.ZipFile(resources_zip, 'r') as zip_ref:
            zip_ref.extractall(resources_dir)
    else:
        print(f"Error: Resource file not found {resources_zip}")

    # Check if gcmcH.py exists
    gcmch_path = os.path.join(os.path.dirname(__file__), 'python', 'gcmcH.py')

    update_values_path = os.path.join(os.path.dirname(__file__), 'scripts', 'update_values.py')
    spec = importlib.util.spec_from_file_location("update_values", update_values_path)
    update_values = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(update_values)
    
    # Run the update_values_py function
    header_file = os.path.join(os.path.dirname(__file__), 'cpp', 'gcmc.h')
    update_values.update_values_py(header_file, gcmch_path)
    print(f"Generated {gcmch_path}")



from .python import GCMC, packages


