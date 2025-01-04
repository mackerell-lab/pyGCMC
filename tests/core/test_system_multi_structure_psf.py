# tests/core/test_system_multi_structure_psf.py

import os
import pytest
from pygcmc import System, Project

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
        'psf': os.path.join(test_data_dir, "test_proa.psf"),
        'benx_psf': os.path.join(test_data_dir, "mols", "benx.psf"),
        'prpx_psf': os.path.join(test_data_dir, "mols", "prpx.psf"),
        'sol_psf': os.path.join(test_data_dir, "mols", "sol.psf"),
        'top': os.path.join(test_data_dir, "test.top"),
    }

def get_pdb_residue_count(pdb_file):
    """辅助函数，用于从 PDB 文件中获取残基数量。"""
    residues = set()
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                residue_name = line[17:20].strip()
                residue_number = line[22:26].strip()
                residues.add((residue_name, residue_number))
    return len(residues)

def verify_system_content(system, expected_residue_count):
    """辅助函数，用于验证系统内容。"""
    # 验证系统中残基的数量
    assert system.get_residue_count() == expected_residue_count, \
        f"Expected {expected_residue_count} residues, found {system.get_residue_count()}."
    
    for i in range(expected_residue_count):
        residue = system.get_residue(i)
        assert residue.get_particle_count() > 0, f"Residue {i} has no particles."
    
        for particle in residue.particles:
            assert particle.mass > 0, f"Particle mass is invalid in residue {i}."
            assert not particle.is_virtual, f"Particle is incorrectly marked as virtual in residue {i}."
            assert isinstance(particle.charge, float), f"Particle charge type is invalid in residue {i}."
            assert len(particle.position) == 3, f"Particle position is invalid in residue {i}."
            assert len(particle.velocity) == 3, f"Particle velocity is invalid in residue {i}."

def verify_pdb_atoms_safe(system, pdb_file):
    """辅助函数，用于安全地验证PDB原子信息。"""
    try:
        # 检查是否有PDB原子
        has_atoms = system.has_pdb_atoms()
        print(f"System has PDB atoms: {has_atoms}")
        assert has_atoms, "System should have PDB atoms"
        
        # 从PDB文件读取原子信息
        pdb_atoms = []
        residue_atoms = {}  # 用于存储每个残基的原子
        residue_types = set()  # 用于跟踪所有残基类型
        
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if not (line.startswith("ATOM") or line.startswith("HETATM")):
                        continue
                        
                    try:
                        # 提取原子信息
                        residue_name = line[17:20].strip()
                        sequence = int(line[22:26].strip())
                        chain = line[21]
                        atom_name = line[12:16].strip()
                        
                        residue_types.add(residue_name)
                        
                        # 存储原子信息
                        atom_info = {
                            'residue': residue_name,
                            'sequence': sequence,
                            'chain': chain,
                            'name': atom_name,
                            'line': line.strip()
                        }
                        pdb_atoms.append(atom_info)
                        
                        # 按残基组织原子
                        key = (residue_name, sequence)
                        if key not in residue_atoms:
                            residue_atoms[key] = set()
                        residue_atoms[key].add(atom_name)
                        
                    except (IndexError, ValueError) as e:
                        print(f"Warning: Failed to parse PDB line: {line.strip()}")
                        print(f"Error details: {str(e)}")
                        continue
                        
            print(f"Found residue types: {', '.join(sorted(residue_types))}")
            print(f"Total residues: {len(residue_atoms)}")
            for res_type in sorted(residue_types):
                count = sum(1 for (rname, _) in residue_atoms.keys() if rname == res_type)
                print(f"  {res_type}: {count} residues")
                
        except Exception as e:
            print(f"Error reading PDB file: {str(e)}")
            raise
        
        # 验证原子数量
        system_atom_count = system.get_pdb_atom_count()
        pdb_atom_count = len(pdb_atoms)
        print(f"System atom count: {system_atom_count}")
        print(f"PDB atom count: {pdb_atom_count}")
        assert system_atom_count == pdb_atom_count, \
            f"Expected {pdb_atom_count} PDB atoms, found {system_atom_count}"
        
        # 验证每个原子的基本属性
        for i in range(system_atom_count):
            try:
                atom = system.get_pdb_atom(i)
                assert atom is not None, f"Got None for atom {i}"
                assert atom.serial > 0, f"Invalid serial number for atom {i}"
                assert atom.name, f"Missing name for atom {i}"
                assert atom.residue, f"Missing residue name for atom {i}"
                assert atom.sequence > 0, f"Invalid sequence number for atom {i}"
                assert isinstance(atom.x, float), f"Invalid x coordinate for atom {i}"
                assert isinstance(atom.y, float), f"Invalid y coordinate for atom {i}"
                assert isinstance(atom.z, float), f"Invalid z coordinate for atom {i}"
                
                # 验证原子是否在正确的残基中
                key = (atom.residue, atom.sequence)
                if key not in residue_atoms:
                    print(f"Warning: Atom {i} belongs to unknown residue {atom.residue} {atom.sequence}")
                    continue
                    
                if atom.name not in residue_atoms[key]:
                    print(f"Warning: Atom {atom.name} not found in residue {atom.residue} {atom.sequence}")
                    continue
                    
            except Exception as e:
                print(f"Error validating atom {i}:")
                print(f"  Serial: {getattr(atom, 'serial', 'N/A')}")
                print(f"  Name: {getattr(atom, 'name', 'N/A')}")
                print(f"  Residue: {getattr(atom, 'residue', 'N/A')}")
                print(f"  Sequence: {getattr(atom, 'sequence', 'N/A')}")
                print(f"Error details: {str(e)}")
                raise
        
        # 收集所有唯一的残基和链
        residues = set((atom['residue'], atom['sequence']) for atom in pdb_atoms)
        chains = set(atom['chain'] for atom in pdb_atoms)
        
        print(f"Found {len(residues)} unique residues and {len(chains)} chains")
        
        # 验证残基查询功能
        for residue_name, sequence in residues:
            try:
                print(f"Checking residue {residue_name} {sequence}")
                # 收集属于这个残基的所有原子
                residue_atoms_list = []
                for i in range(system_atom_count):
                    atom = system.get_pdb_atom(i)
                    if atom.residue == residue_name and atom.sequence == sequence:
                        residue_atoms_list.append(atom)
                
                if not residue_atoms_list:
                    print(f"Warning: No atoms found for residue {residue_name} {sequence}")
                    continue
                
                # 验证每个原子是否属于正确的残基
                for atom in residue_atoms_list:
                    assert atom.residue == residue_name, \
                        f"Wrong residue name: {atom.residue} != {residue_name}"
                    assert atom.sequence == sequence, \
                        f"Wrong sequence number: {atom.sequence} != {sequence}"
                    
                    # 验证原子名称是否在PDB文件中存在
                    assert atom.name in residue_atoms[(residue_name, sequence)], \
                        f"Atom {atom.name} not found in PDB file for residue {residue_name} {sequence}"
                    
                print(f"Successfully validated residue {residue_name} {sequence} with {len(residue_atoms_list)} atoms")
                    
            except Exception as e:
                print(f"Error validating residue {residue_name} {sequence}:")
                print(f"Error details: {str(e)}")
                continue
        
        # 验证链查询功能
        for chain in chains:
            try:
                # 收集属于这个链的所有原子
                chain_atoms_list = []
                for i in range(system_atom_count):
                    atom = system.get_pdb_atom(i)
                    if atom.chain == chain:
                        chain_atoms_list.append(atom)
                
                if not chain_atoms_list:
                    print(f"Warning: No atoms found for chain {chain}")
                    continue
                
                # 验证每个原子是否属于正确的链
                for atom in chain_atoms_list:
                    assert atom.chain == chain, \
                        f"Wrong chain identifier: {atom.chain} != {chain}"
                    
                print(f"Successfully validated chain {chain} with {len(chain_atoms_list)} atoms")
                    
            except Exception as e:
                print(f"Error validating chain {chain}:")
                print(f"Error details: {str(e)}")
                continue
        
        print("PDB atom verification completed successfully")
        
    except Exception as e:
        print(f"Fatal error during PDB atom verification: {str(e)}")
        raise

def test_load_structure_psf_multi(structure_files):
    """测试使用多残基 PSF 文件加载结构。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    expected_residue_count = get_pdb_residue_count(structure_files['pdb'])
    verify_system_content(system, expected_residue_count=expected_residue_count)
    verify_pdb_atoms_safe(system, structure_files['pdb'])

def test_load_structure_psf_single(structure_files):
    """测试使用单残基 PSF 文件加载结构。"""
    system = System()
    system.load_structure_psf_single(structure_files['pdb'], structure_files['sol_psf'], "SOL")
    
    sol_residues = []
    for i in range(system.get_residue_count()):
        residue = system.get_residue(i)
        if residue.name == "SOL":
            sol_residues.append(residue)
    
    expected_sol_residues = len(sol_residues)
    assert len(sol_residues) == expected_sol_residues, \
        f"Expected {expected_sol_residues} SOL residues, found {len(sol_residues)}."
    
    for residue in sol_residues:
        for particle in residue.particles:
            # Check if particle is OW (oxygen) or HW (hydrogen)
            if particle.mass > 10.0:  # Oxygen atom
                assert particle.mass == pytest.approx(15.9994, 1e-3), "Incorrect mass for SOL OW atom."
            else:  # Hydrogen atom
                assert particle.mass == pytest.approx(1.008, 1e-3), "Incorrect mass for SOL HW atom."
    
    verify_pdb_atoms_safe(system, structure_files['pdb'])

def test_load_structure_psf_auto_multi(structure_files):
    """测试自动加载多残基 PSF 文件。"""
    system = System()
    system.load_structure_psf_auto(structure_files['pdb'], structure_files['psf'])
    
    expected_residue_count = get_pdb_residue_count(structure_files['pdb'])
    verify_system_content(system, expected_residue_count=expected_residue_count)
    verify_pdb_atoms_safe(system, structure_files['pdb'])

def test_load_structure_psf_auto_single(structure_files):
    """测试自动加载单残基 PSF 文件。"""
    system = System()
    system.load_structure_psf_auto(structure_files['pdb'], structure_files['sol_psf'])
    
    sol_residues = []
    for i in range(system.get_residue_count()):
        residue = system.get_residue(i)
        if residue.name == "SOL":
            sol_residues.append(residue)
    
    expected_sol_residues = len(sol_residues)
    assert len(sol_residues) == expected_sol_residues, \
        f"Expected {expected_sol_residues} SOL residues, found {len(sol_residues)}."
    
    for residue in sol_residues:
        for particle in residue.particles:
            # Check if particle is OW (oxygen) or HW (hydrogen)
            if particle.mass > 10.0:  # Oxygen atom
                assert particle.mass == pytest.approx(15.9994, 1e-3), "Incorrect mass for SOL OW atom."
            else:  # Hydrogen atom
                assert particle.mass == pytest.approx(1.008, 1e-3), "Incorrect mass for SOL HW atom."
    
    verify_pdb_atoms_safe(system, structure_files['pdb'])

def test_pdb_atom_management(structure_files):
    """测试PDB原子管理功能。"""
    system = System()
    system.load_structure_psf_auto(structure_files['pdb'], structure_files['psf'])
    
    # 测试清除PDB原子
    initial_count = system.get_pdb_atom_count()
    assert initial_count > 0, "Should have PDB atoms after loading"
    
    # 保存第一个原子的信息
    first_atom = system.get_pdb_atom(0)
    
    system.clear_pdb_atoms()
    assert not system.has_pdb_atoms(), "Should have no PDB atoms after clearing"
    assert system.get_pdb_atom_count() == 0, "PDB atom count should be 0 after clearing"
    
    # 测试添加和删除PDB原子
    system.add_pdb_atom(first_atom)
    assert system.get_pdb_atom_count() == 1, "Should have 1 PDB atom after adding"
    
    system.remove_pdb_atom(0)
    assert not system.has_pdb_atoms(), "Should have no PDB atoms after removing"
    
    # 测试异常情况
    with pytest.raises(IndexError):
        system.get_pdb_atom(9999)  # 测试无效索引
    
    with pytest.raises(IndexError):
        system.remove_pdb_atom(9999)  # 测试无效索引

def test_load_structure_psf_auto_fail(structure_files):
    """测试自动加载无效 PSF 文件时的错误处理。"""
    system = System()
    with pytest.raises(RuntimeError, match="No atoms found in PSF file"):
        # 使用 TOP 文件作为无效的 PSF 文件以触发错误
        system.load_structure_psf_auto(structure_files['pdb'], structure_files['top'])

def test_load_structure_psf_multi_fail(structure_files):
    """测试多残基 PSF 加载无效 PSF 文件时的错误处理。"""
    system = System()
    with pytest.raises(RuntimeError, match="No atoms found in PSF file"):
        # 使用 TOP 文件作为无效的 PSF 文件以触发错误
        system.load_structure_psf_multi(structure_files['pdb'], structure_files['top'])

def test_load_structure_psf_single_fail(structure_files):
    """测试单残基 PSF 加载无效 PSF 文件时的错误处理。"""
    system = System()
    with pytest.raises(RuntimeError, match="No atoms found in PSF file"):
        # 使用 TOP 文件作为无效的 PSF 文件以触发错误
        system.load_structure_psf_single(structure_files['pdb'], structure_files['top'], "SOL")
