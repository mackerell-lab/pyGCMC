import re

class Atom:
    def __init__(self, serial, name, residue, sequence, x, y, z, s):
        self.serial = serial
        self.name = name
        self.residue = residue
        self.sequence = sequence
        self.x = x
        self.y = y
        self.z = z
        self.s = s


def read_pdb(pdb_file):
    '''Read a PDB file and return the crystal xyz and a list of atoms'''

    s = open(pdb_file, 'r').read()
    p = s.strip().split('\n')
    cryst = None
    for line in p:
        if line.lower().startswith('cryst1 '):
            cryst = line.split()[1:4]
            cryst = [float(i) for i in cryst]
            break

    atoms = []
    for line in p:
        if (line.lower().startswith('atom ') or line.lower().startswith('hetatm')) and len(line) > 54:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            serial = line[6:11].strip()
            name = line[12:16].strip()
            residue = line[17:20].strip()
            sequence = line[22:26].strip()
            atom = Atom(serial, name, residue, sequence, x, y, z, line)
            atoms.append(atom)
    return cryst,atoms



def read_itp(itp_file):
    '''Read a Gromacs .itp file and return a list of atoms'''

    s = open(itp_file, 'r').read()

    blocks = {}
    current_block = "title"
    blocks[current_block] = []
    for line in s.split('\n'):
        line = re.sub(r';.*', '', line)
        line = re.sub(r'#.*', '', line)
        match = re.match(r'\[(.*)\]', line)
        if match:
            current_block = match.group(1).strip()
            if current_block not in blocks:
                blocks[current_block] = []
        else:
            if current_block and line.strip():
                blocks[current_block].append(line.strip())


    return blocks

def read_top(top_file):
    '''Read a Gromacs .top file and return a list of atoms'''

    s = open(top_file, 'r').read()


    blocks = {}
    current_block = "title"
    blocks[current_block] = []
    for line in s.split('\n'):
        line = re.sub(r';.*', '', line)
        match = re.match(r'\[(.*)\]', line)
        if match:
            current_block = match.group(1).strip()
            if current_block not in blocks:
                blocks[current_block] = []
        else:
            if current_block and line.strip():
                blocks[current_block].append(line.strip())
    
    content = re.sub(r';.*', '', s)
    content = re.sub(r'#ifdef .*?#endif', '', content, flags=re.DOTALL)

    pattern = r'^\s*#include "(.*?)"\s*$'
    matches = re.findall(pattern, content, re.MULTILINE)
    matches = [i.strip() for i in matches]

    itps = []
    for file in matches:
        if file.endswith('.itp'):
            itps += [read_itp(file)]
            

    print(blocks.keys())
    print(blocks['moleculetype'])
    print(matches)
    print(itps)

