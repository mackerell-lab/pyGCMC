import re

class Atom:
    def __init__(self, serial, name, residue, sequence, x, y, z):
        self.serial = serial
        self.name = name
        self.residue = residue
        self.sequence = sequence
        self.x = x
        self.y = y
        self.z = z


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
            atom = Atom(serial, name, residue, sequence, x, y, z)
            atoms.append(atom)
    return cryst,atoms



def read_itp(itp_file):
    '''Read a Gromacs .itp file and return a list of atoms'''

    s = open(itp_file, 'r').read()

    # blocks = {}
    # current_block = "title"
    # blocks[current_block] = []
    # for line in s.split('\n'):
    #     line = re.sub(r';.*', '', line)
    #     line = re.sub(r'#.*', '', line)
    #     match = re.match(r'\[(.*)\]', line)
    #     if match:
    #         current_block = match.group(1).strip()
    #         if current_block not in blocks:
    #             blocks[current_block] = []
    #     else:
    #         if current_block and line.strip():
    #             blocks[current_block].append(line.strip())
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

def read_top(top_file):
    '''Read a Gromacs .top file and return a list of atoms'''

    s = open(top_file, 'r').read()


    # blocks = {}
    # current_block = "title"
    # blocks[current_block] = []
    # for line in s.split('\n'):
    #     line = re.sub(r';.*', '', line)
    #     match = re.match(r'\[(.*)\]', line)
    #     if match:
    #         current_block = match.group(1).strip()
    #         if current_block not in blocks:
    #             blocks[current_block] = []
    #     else:
    #         if current_block and line.strip():
    #             blocks[current_block].append(line.strip())
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

    itps = []
    for file in matches:
        if file.endswith('.itp'):
            itps += [read_itp(file)]
    
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

    for blocknum in range(len(blocks)):
        block = blocks[blocknum]
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
                            # else:
                            #     print(blocknum,line)
                            #     break
                    break

    for itp in itps:
        for blocknum in range(len(itp)):
            block = itp[blocknum]
            if block[0] == 'moleculetype':
                for line in block[1:]:
                    line = line.strip()
                    if line:
                        pattern = r'(\w+)\s+(\d+)'
                        match = re.match(pattern, line)
                        if match:
                            moleculetypes.append([match.group(1)])
                            break
                for atomblock in itp[blocknum+1:]:
                    if atomblock[0] == 'atoms':
                        for line in atomblock[1:]:
                            line = line.strip()
                            if line:
                                pattern = r'(\d+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)\s+(\d+)\s+([0-9\.\-]+)'
                                match = re.match(pattern, line)
                                if match:
                                    moleculetypes[-1].append([match.group(2),match.group(3),match.group(4),match.group(5),match.group(7)])
                                # else:
                                #     print(atomblock[0],line)
                                #     break
                        break
    
    atoms_top = []
    for molecule in molecules:
        for moleculetype in moleculetypes:
            if molecule[0] == moleculetype[0]:
                # print(molecule[0])
                # print(len(moleculetype))
                for i in range(int(molecule[1])):
                    atoms_top+=moleculetype[1:]
                break
    # print(moleculetype)
    # print(atoms_top[233336])
    return atoms_top
    # for i in moleculetypes:
    #     print(i[0])
    
    # print(molecules)
    # print([i[0] for i in blocks])      
    # # print(blocks['molecules'])
    # # print(blocks['moleculetype'])
    # for block in itps:
    #     print([i[0] for i in block])    
    #     for i in block:
    #         if i[0] == 'moleculetype':
    #             print(i)

