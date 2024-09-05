"""
    © Copyright 2023 - University of Maryland, Baltimore   All Rights Reserved    
    	Mingtian Zhao, Alexander D. MacKerell Jr.        
    E-mail: 
    	zhaomt@outerbanks.umaryland.edu
    	alex@outerbanks.umaryland.edu
"""

import re

# Atom class to represent individual atoms
class Atom:
    def __init__(self, serial, name, residue, sequence, x, y, z):
        self.serial = serial
        self.name = name
        self.residue = residue
        self.sequence = sequence
        self.x = x
        self.y = y
        self.z = z
        self.s = None

# Read PDB file and extract crystal information and atom data
def read_pdb(pdb_file):
    with open(pdb_file, 'r') as f:
        s = f.read()
    p = s.strip().split('\n')
    
    # Extract crystal information
    cryst = None
    for line in p:
        if line.lower().startswith('cryst1 '):
            cryst = [float(i) for i in line.split()[1:4]]
            break

    # Extract atom information
    atoms = []
    for line in p:
        if (line.lower().startswith('atom ') or line.lower().startswith('hetatm')) and len(line) > 54:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            serial = line[6:11].strip()
            name = line[12:16].strip()
            residue = line[17:20].strip()
            sequence = line[22:26].strip()
            atom = Atom(serial, name, residue, sequence, x, y, z)
            atom.s = line
            atoms.append(atom)

    # Calculate crystal dimensions if not provided
    if cryst is None:
        x = [atom.x for atom in atoms]
        y = [atom.y for atom in atoms]
        z = [atom.z for atom in atoms]
        cryst = [max(x)-min(x), max(y)-min(y), max(z)-min(z)]
    
    return cryst, atoms, s

# Read and parse ITP file blocks
def read_itp_block(itp_file):
    with open(itp_file, 'r') as f:
        s = f.read()

    blocks = []
    current_block = "title"
    blocks.append([current_block])
    for line in s.split('\n'):
        line = re.sub(r';.*', '', line)
        line = re.sub(r'#.*', '', line)
        match = re.match(r'\[(.*)\]', line)
        if match:
            current_block = match.group(1).strip()
            blocks.append([current_block])
        else:
            if current_block and line.strip():
                blocks[-1].append(line.strip())

    return blocks

# Read ITP file and extract molecule types and atom information
def read_itp(itp_file):
    blocks = read_itp_block(itp_file)
    moleculetypes = []
    for blocknum, block in enumerate(blocks):
        if block[0] == 'moleculetype':
            for line in block[1:]:
                line = line.strip()
                if line:
                    pattern = r'(\w+)\s+(\d+)'
                    match = re.match(pattern, line)
                    if match:
                        moleculetypes.append([match.group(1)])
                        break
            for atomblock in blocks[blocknum+1:]:
                if atomblock[0] == 'atoms':
                    for line in atomblock[1:]:
                        line = line.strip()
                        if line:
                            pattern = r'(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+([0-9\.\-]+)'
                            match = re.match(pattern, line)
                            if match:
                                moleculetypes[-1].append([match.group(2),match.group(3),match.group(4),match.group(5),match.group(7)])
                    break
    return moleculetypes[-1][1:]

# Read PSF file and extract atom information
def read_psf(psf_file):
    with open(psf_file, 'r') as f:
        s = f.read()

    p = s.strip().split('\n')[1:]  # skip the first line
    p = [line.split('*', 1)[0].split('!', 1)[0].split('#', 1)[0] for line in p]

    for i, line in enumerate(p):
        try:
            lineN = int(line)
            linei = i
            break
        except:
            pass

    p = p[linei + lineN + 1:]

    for i, line in enumerate(p):
        try:
            lineN = int(line)
            linei = i
            break
        except:
            pass

    p = p[linei + 1:]

    atoms_top = []
    n = 0
    for line in p:
        line = line.strip()
        pattern = r"^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)"
        match = re.match(pattern, line)
        if match:
            n += 1
            if n > lineN:
                break
            atoms_top.append([match.group(6),match.group(3),match.group(4),match.group(5),match.group(7)])
    
    return atoms_top

# Read TOP file and extract molecule and atom information
def read_top(top_file):
    with open(top_file, 'r') as f:
        s = f.read()

    blocks = []
    current_block = "title"
    blocks.append([current_block])
    for line in s.split('\n'):
        line = re.sub(r';.*', '', line)
        line = re.sub(r'#.*', '', line)
        match = re.match(r'\[(.*)\]', line)
        if match:
            current_block = match.group(1).strip()
            blocks.append([current_block])
        else:
            if current_block and line.strip():
                blocks[-1].append(line.strip())
    
    content = re.sub(r';.*', '', s)
    content = re.sub(r'#ifdef .*?#endif', '', content, flags=re.DOTALL)

    pattern = r'^\s*#include "(.*?)"\s*$'
    matches = re.findall(pattern, content, re.MULTILINE)
    matches = [i.strip() for i in matches]

    itps = [read_itp_block(file) for file in matches if file.endswith('.itp')]
    
    molecules = []
    for block in blocks:
        if block[0] == 'molecules':
            for line in block[1:]:
                line = line.strip()
                if line:
                    pattern = r'(\w+)\s+(\d+)'
                    match = re.match(pattern, line)
                    if match:
                        molecules.append([match.group(1),match.group(2)])
            break
    
    moleculetypes = []

    for blocknum, block in enumerate(blocks):
        if block[0] == 'moleculetype':
            moleculetypes.append(process_moleculetype(block, blocks[blocknum+1:]))

    for itp in itps:
        for blocknum, block in enumerate(itp):
            if block[0] == 'moleculetype':
                moleculetypes.append(process_moleculetype(block, itp[blocknum+1:]))
    
    atoms_top = []
    for molecule in molecules:
        for moleculetype in moleculetypes:
            if molecule[0] == moleculetype[0]:
                atoms_top += moleculetype[1:] * int(molecule[1])
                break

    return atoms_top, s

# Helper function to process moleculetype blocks
def process_moleculetype(block, following_blocks):
    moleculetype = []
    for line in block[1:]:
        line = line.strip()
        if line:
            pattern = r'(\w+)\s+(\d+)'
            match = re.match(pattern, line)
            if match:
                moleculetype.append(match.group(1))
                break
    for atomblock in following_blocks:
        if atomblock[0] == 'atoms':
            for line in atomblock[1:]:
                line = line.strip()
                if line:
                    pattern = r'(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+([0-9\.\-]+)'
                    match = re.match(pattern, line)
                    if match:
                        moleculetype.append([match.group(2),match.group(3),match.group(4),match.group(5),match.group(7)])
            break
    return moleculetype

# Read FF file and extract non-bonded parameters
def read_ff(ff_file):
    nb_dict = {}
    nbfix_dict = {}
    blocks = read_itp_block(ff_file)

    for block in blocks:
        if block[0] == 'atomtypes':
            for line in block[1:]:
                line = line.strip()
                if line:
                    pattern = r'(\S+)\s+(\S+)\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s+(\S+)\s+([0-9\.\-]+)+\s+([0-9\.\-]+)'
                    match = re.match(pattern, line)
                    if match:
                        nb_dict[match.group(1)] = [float(match.group(6)),float(match.group(7))]
    
        if block[0] in ['nonbond_params', 'pairtypes']:
            for line in block[1:]:
                line = line.strip()
                if line:
                    pattern = r'(\S+)\s+(\S+)\s+([0-9\.\-]+)\s+([0-9\.\-]+)\s+([0-9\.\-]+)'
                    match = re.match(pattern, line)
                    if match:
                        nbfix_dict[(match.group(1),match.group(2))] = [float(match.group(4)),float(match.group(5))]

    return nb_dict, nbfix_dict
