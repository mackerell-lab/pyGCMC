# tests/core/test_system_pdbatom.py

import os
import pytest
import math
from pygcmc import System

def calculate_distance(atom1, atom2):
    """计算两个原子之间的距离。"""
    return math.sqrt(
        (atom1.x - atom2.x)**2 +
        (atom1.y - atom2.y)**2 +
        (atom1.z - atom2.z)**2
    )

def get_atom_by_name(atoms, name):
    """从原子列表中获取指定名称的原子。"""
    for atom in atoms:
        if atom.name == name:
            return atom
    return None

def is_standard_amino_acid(residue_name):
    """检查是否为标准氨基酸。"""
    return residue_name in {
        "ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "TYR", "HIS",
        "LYS", "ARG", "GLU", "ASP", "GLN", "ASN", "MET", "PRO",
        "THR", "SER", "CYS", "GLY"
    }

def is_solvent(residue_name):
    """检查是否为溶剂分子。"""
    return residue_name in {"SOL", "WAT", "HOH", "TIP3"}

@pytest.fixture
def test_data_dir():
    """获取测试数据目录的路径。"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

@pytest.fixture
def structure_files(test_data_dir):
    """获取结构文件。"""
    return {
        'pdb': os.path.join(test_data_dir, "test.pdb"),
        'psf': os.path.join(test_data_dir, "test_proa.psf")
    }

def test_system_initialization(structure_files):
    """测试系统初始化时PDB和PSF数据的加载。"""
    system = System(pdb=structure_files['pdb'], psf=structure_files['psf'])
    
    # 验证系统成功加载了原子
    assert system.has_pdb_atoms(), "System should have PDB atoms after initialization"
    assert system.get_pdb_atom_count() > 0, "System should have non-zero PDB atoms"
    
    # 验证第一个原子的完整性
    first_atom = system.get_pdb_atom(0)
    assert first_atom is not None, "First atom should not be None"
    assert first_atom.serial == 1, "First atom should have serial number 1"
    assert first_atom.name == "N", "First atom should be N"
    assert first_atom.residue == "ALA", "First atom should be in ALA"
    assert first_atom.sequence == 7, "First atom should be in residue 7"

def test_pdb_data_integrity(structure_files):
    """测试PDB数据的完整性。"""
    system = System(pdb=structure_files['pdb'], psf=structure_files['psf'])
    
    # 读取原始PDB文件以获取参考数据
    pdb_atoms = []
    with open(structure_files['pdb'], 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                # 解析PDB行
                serial = int(line[6:11])
                name = line[12:16].strip()
                residue = line[17:20].strip()
                sequence = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                
                # 只存储蛋白质原子的数据
                if is_standard_amino_acid(residue):
                    pdb_atoms.append({
                        'serial': serial,
                        'name': name,
                        'residue': residue,
                        'sequence': sequence,
                        'x': x,
                        'y': y,
                        'z': z
                    })
    
    # 验证蛋白质原子的属性
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        if is_standard_amino_acid(atom.residue):
            # 在参考数据中查找匹配的原子
            ref_atom = next((a for a in pdb_atoms 
                            if a['serial'] == atom.serial), None)
            assert ref_atom is not None, f"No reference atom found for serial {atom.serial}"
            
            # 验证原子属性
            assert atom.name == ref_atom['name'], \
                f"Wrong name for atom {i}: {atom.name} != {ref_atom['name']}"
            assert atom.residue == ref_atom['residue'], \
                f"Wrong residue for atom {i}: {atom.residue} != {ref_atom['residue']}"
            assert atom.sequence == ref_atom['sequence'], \
                f"Wrong sequence for atom {i}: {atom.sequence} != {ref_atom['sequence']}"
            assert atom.x == pytest.approx(ref_atom['x']), \
                f"Wrong x coordinate for atom {i}"
            assert atom.y == pytest.approx(ref_atom['y']), \
                f"Wrong y coordinate for atom {i}"
            assert atom.z == pytest.approx(ref_atom['z']), \
                f"Wrong z coordinate for atom {i}"

def test_psf_data_integrity(structure_files):
    """测试PSF数据的完整性。"""
    system = System(pdb=structure_files['pdb'], psf=structure_files['psf'])
    
    # 读取原始PSF文件以获取参考数据
    psf_atoms = []
    reading_atoms = False
    with open(structure_files['psf'], 'r') as f:
        for line in f:
            if "!NATOM" in line:
                reading_atoms = True
                continue
            if reading_atoms:
                if line.strip() == "":
                    break
                parts = line.split()
                if len(parts) >= 7:
                    psf_atoms.append({
                        'id': int(parts[0]),
                        'segment': parts[1],
                        'residue_id': int(parts[2]),
                        'residue': parts[3],
                        'name': parts[4],
                        'type': parts[5],
                        'charge': float(parts[6]),
                        'mass': float(parts[7])
                    })
    
    # 验证拓扑信息（只验证蛋白质原子）
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        if is_standard_amino_acid(atom.residue):
            # 查找对应的PSF原子（考虑残基ID和原子名称）
            psf_atom = next((a for a in psf_atoms 
                            if a['name'] == atom.name and 
                            a['residue'] == atom.residue and
                            a['residue_id'] == atom.sequence), None)
            if psf_atom:
                if atom.topo_type:  # 只在有拓扑类型时验证
                    assert atom.topo_type == psf_atom['type'], \
                        f"Wrong topology type for atom {i}"
                if not math.isnan(atom.topo_charge):  # 只在有拓扑电荷时验证
                    assert atom.topo_charge == pytest.approx(psf_atom['charge']), \
                        f"Wrong topology charge for atom {i}"
                if not math.isnan(atom.topo_mass):  # 只在有拓扑质量时验证
                    assert atom.topo_mass == pytest.approx(psf_atom['mass']), \
                        f"Wrong topology mass for atom {i}"

def test_residue_consistency(structure_files):
    """测试残基信息的一致性。"""
    system = System(pdb=structure_files['pdb'], psf=structure_files['psf'])
    
    # 收集所有唯一的残基
    residues = set()
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        residues.add((atom.residue, atom.sequence))
    
    # 验证每个残基的完整性
    for residue_name, sequence in sorted(residues):
        atoms = system.get_pdb_atoms_by_residue_sequence(residue_name, sequence)
        assert len(atoms) > 0, f"No atoms found for residue {residue_name} {sequence}"
        
        # 验证所有原子属于同一个残基
        for atom in atoms:
            assert atom.residue == residue_name, \
                f"Wrong residue name: {atom.residue} != {residue_name}"
            assert atom.sequence == sequence, \
                f"Wrong sequence number: {atom.sequence} != {sequence}"
        
        # 验证主链原子的存在性（对于标准氨基酸，排除末端残基和特殊残基）
        if (residue_name in {"ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "TYR", "HIS",
                          "LYS", "ARG", "GLU", "ASP", "GLN", "ASN", "MET", "PRO",
                          "THR", "SER", "CYS", "GLY"} and
            sequence not in {7, max(s for _, s in residues)}):  # 排除首尾残基
            backbone_atoms = {"N", "CA", "C"}  # 只检查必需的主链原子
            atom_names = {atom.name for atom in atoms}
            for backbone_atom in backbone_atoms:
                assert backbone_atom in atom_names, \
                    f"Missing backbone atom {backbone_atom} in {residue_name} {sequence}"

def test_chain_consistency(structure_files):
    """测试链信息的一致性。"""
    system = System(pdb=structure_files['pdb'], psf=structure_files['psf'])
    
    # 收集所有链
    chains = set()
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        chains.add(atom.chain)
    
    # 验证每条链
    for chain in chains:
        atoms = system.get_pdb_atoms_by_chain(chain)
        assert len(atoms) > 0, f"No atoms found for chain {chain}"
        
        # 验证所有原子属于同一条链
        for atom in atoms:
            assert atom.chain == chain, \
                f"Wrong chain identifier: {atom.chain} != {chain}"
        
        # 验证残基序号的连续性
        sequences = sorted(list(set(atom.sequence for atom in atoms)))
        for i in range(len(sequences) - 1):
            assert sequences[i+1] > sequences[i], \
                f"Non-sequential residue numbers in chain {chain}"

def test_topology_completeness(structure_files):
    """测试拓扑信息的完整性。"""
    system = System(pdb=structure_files['pdb'], psf=structure_files['psf'])
    
    # 检查所有原子是否都有拓扑信息
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        
        # 验证拓扑属性的存在性
        assert hasattr(atom, 'topo_type'), f"No topology type attribute for atom {i}"
        assert hasattr(atom, 'topo_charge'), f"No topology charge attribute for atom {i}"
        assert hasattr(atom, 'topo_mass'), f"No topology mass attribute for atom {i}"
        
        # 如果有拓扑信息，验证其合理性
        if atom.topo_type:
            assert len(atom.topo_type) > 0, f"Empty topology type for atom {i}"
        if not math.isnan(atom.topo_charge):
            assert -2.0 <= atom.topo_charge <= 2.0, \
                f"Unreasonable topology charge for atom {i}: {atom.topo_charge}"
        if not math.isnan(atom.topo_mass):
            assert 1.0 <= atom.topo_mass <= 200.0, \
                f"Unreasonable topology mass for atom {i}: {atom.topo_mass}"

def test_coordinate_validity(structure_files):
    """测试坐标的有效性。"""
    system = System(pdb=structure_files['pdb'], psf=structure_files['psf'])
    
    # 检查所有原子的坐标
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        
        # 验证坐标是有限数
        assert all(map(math.isfinite, [atom.x, atom.y, atom.z])), \
            f"Invalid coordinates for atom {i}"
        
        # 验证坐标在合理范围内（根据CRYST1记录）
        assert 0 <= atom.x <= 127.022, f"X coordinate out of box for atom {i}"
        assert 0 <= atom.y <= 133.419, f"Y coordinate out of box for atom {i}"
        assert 0 <= atom.z <= 132.854, f"Z coordinate out of box for atom {i}"

def test_special_positions(structure_files):
    """测试特殊位置原子的处理。"""
    system = System(pdb=structure_files['pdb'], psf=structure_files['psf'])
    
    # 获取第一个和最后一个原子
    first_atom = system.get_pdb_atom(0)
    last_atom = system.get_pdb_atom(system.get_pdb_atom_count() - 1)
    
    # 验证第一个原子（通常是N端）
    assert first_atom.name == "N", "First atom should be N"
    if first_atom.topo_type:  # 只在有拓扑类型时验证
        assert first_atom.topo_type == "NH3", "First atom should have NH3 type"
        assert first_atom.topo_charge == pytest.approx(-0.300000), \
            "First atom should have correct charge"
    
    # 验证最后一个原子的有效性
    assert last_atom.serial > 0, "Last atom should have valid serial number"
    assert last_atom.name, "Last atom should have valid name"
    assert last_atom.residue, "Last atom should have valid residue name"
    assert last_atom.sequence > 0, "Last atom should have valid sequence number"
    
    # 验证N端残基
    n_term_atoms = system.get_pdb_atoms_by_residue_sequence("ALA", 7)  # N端
    assert len(n_term_atoms) > 0, "N-terminal residue should have atoms"
    n_term_atom_names = {atom.name for atom in n_term_atoms}
    assert "N" in n_term_atom_names, "N-terminal residue should have N"
    assert any(name.startswith("H") for name in n_term_atom_names), \
        "N-terminal residue should have hydrogen atoms"
    
    # 获取最后一个蛋白质残基（排除溶剂分子和小分子）
    protein_residues = [(atom.residue, atom.sequence) for atom in 
                       [system.get_pdb_atom(i) 
                        for i in range(system.get_pdb_atom_count())]
                       if is_standard_amino_acid(atom.residue)]  # 只考虑标准氨基酸
    if protein_residues:
        last_protein_residue = max(protein_residues, key=lambda x: x[1])
        c_term_atoms = system.get_pdb_atoms_by_residue_sequence(*last_protein_residue)
        c_term_atom_names = {atom.name for atom in c_term_atoms}
        # 只检查C原子，因为O原子可能已经被修饰
        assert "C" in c_term_atom_names, "C-terminal residue should have C"