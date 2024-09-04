"""

    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu

"""


# import numpy as np
# from . import gpu
from .gcmc import *
from .main import main
from .mainOld import main as mainOld

# Load the CHARMM force field parameters
import os
import zipfile
import pkg_resources


# Update the path to toppar_dir
toppar_dir = os.path.join(os.path.dirname(__file__), 'resources', 'toppar')

if not os.path.exists(toppar_dir) or not os.path.exists(os.path.join(toppar_dir, 'silcs.prm')) or os.path.getsize(os.path.join(toppar_dir, 'silcs.prm')) < 10:
    print("CHARMM force field parameters not found.")
    print("Extracting CHARMM force field parameters...")
    # Update the path to resources.zip
    resources_zip = pkg_resources.resource_filename('gcmc', os.path.join('resources', 'resources.zip'))
    with zipfile.ZipFile(resources_zip, 'r') as zip_ref:
        # Update the extraction destination
        zip_ref.extractall(os.path.join(os.path.dirname(__file__), 'resources'))
