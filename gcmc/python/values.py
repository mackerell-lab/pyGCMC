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

# Data type for representing an atom
Atom_dtype = np.dtype([
    ('position', np.float32, (3,)),
    ('charge', np.float32),
    ('type', np.int32)
])

# Data type for representing a residue
Residue_dtype = np.dtype([
    ('position', np.float32, (3,)),
    ('atomNum', np.int32),
    ('atomStart', np.int32)
])

# Data type for representing an array of atoms
AtomArray_dtype = np.dtype([
    ('name', 'S4'),
    ('startRes', np.int32),
    ('muex', np.float32),
    ('conc', np.float32),
    ('confBias', np.int32),
    ('mcTime', np.float32),
    ('totalNum', np.int32),
    ('maxNum', np.int32),
    ('num_atoms', np.int32),
    ('atoms', Atom_dtype, (20,))
])

# Data type for storing simulation information
Info_dtype = np.dtype([
    ('mcsteps', np.int32),
    ('cutoff', np.float32),
    ('grid_dx', np.float32),
    ('startxyz', np.float32, (3,)),
    ('cryst', np.float32, (3,)),
    ('showInfo', np.int32),
    ('cavityFactor', np.float32),
    ('fragTypeNum', np.int32),
    ('totalGridNum', np.int32),
    ('totalResNum', np.int32),
    ('totalAtomNum', np.int32),
    ('ffXNum', np.int32),
    ('ffYNum', np.int32),
    ('PME', np.int32),
    ('seed', np.uint32)
])
