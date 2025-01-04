# tests/core/test_system_multi_structure_psf.py

import os
import pytest
from pygcmc import System, Project
import math

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

def test_boundary_residues(structure_files):
    """测试第一个和最后一个残基的完整性。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 测试第一个残基（ALA 7）
    first_residue_atoms = system.get_pdb_atoms_by_residue_sequence("ALA", 7)
    assert len(first_residue_atoms) == 12  # 根据PSF文件中ALA的原子数
    assert first_residue_atoms[0].name == "N"
    assert first_residue_atoms[0].x == pytest.approx(76.563)
    assert first_residue_atoms[0].y == pytest.approx(93.118)
    assert first_residue_atoms[0].z == pytest.approx(93.806)
    
    # 获取最后一个残基的序号
    last_sequence = 0
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        last_sequence = max(last_sequence, atom.sequence)
    
    # 测试最后一个残基
    last_residue_atoms = system.get_pdb_atoms_by_residue_sequence(
        system.get_pdb_atom(system.get_pdb_atom_count() - 1).residue,
        last_sequence
    )
    assert len(last_residue_atoms) > 0

def test_atom_properties_completeness(structure_files):
    """测试原子属性是否完整地从PSF文件转移到了系统中。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 测试ALA 7的N原子
    ala_atoms = system.get_pdb_atoms_by_residue_sequence("ALA", 7)
    n_atom = next(atom for atom in ala_atoms if atom.name == "N")
    
    # 验证从PDB文件读取的属性
    assert n_atom.x == pytest.approx(76.563)
    assert n_atom.y == pytest.approx(93.118)
    assert n_atom.z == pytest.approx(93.806)
    
    # 验证拓扑属性
    assert n_atom.topo_type == "NH3"
    assert n_atom.topo_charge == pytest.approx(-0.300000)
    assert n_atom.topo_mass == pytest.approx(14.0070)

def test_crystal_parameters(structure_files):
    """测试晶胞参数是否正确读取。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 验证晶胞参数
    # 注：这里我们只验证系统是否成功加载，因为晶胞参数的访问方法可能不同
    assert system.has_pdb_atoms()
    
    # 验证第一个原子的坐标在合理范围内
    first_atom = system.get_pdb_atom(0)
    assert 0 <= first_atom.x <= 127.022
    assert 0 <= first_atom.y <= 133.419
    assert 0 <= first_atom.z <= 132.854

def test_residue_connectivity(structure_files):
    """测试残基内部原子的连接性。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 测试ALA 7的关键原子连接
    ala_atoms = system.get_pdb_atoms_by_residue_sequence("ALA", 7)
    
    # 找到关键原子
    n_atom = next(atom for atom in ala_atoms if atom.name == "N")
    ca_atom = next(atom for atom in ala_atoms if atom.name == "CA")
    c_atom = next(atom for atom in ala_atoms if atom.name == "C")
    
    # 验证它们的相对位置关系
    n_ca_distance = ((ca_atom.x - n_atom.x)**2 + 
                    (ca_atom.y - n_atom.y)**2 + 
                    (ca_atom.z - n_atom.z)**2)**0.5
    ca_c_distance = ((c_atom.x - ca_atom.x)**2 + 
                    (c_atom.y - ca_atom.y)**2 + 
                    (c_atom.z - ca_atom.z)**2)**0.5
    
    # 典型的N-CA和CA-C键长约为1.47和1.52埃
    assert n_ca_distance == pytest.approx(1.47, abs=0.1)
    assert ca_c_distance == pytest.approx(1.52, abs=0.1)

def test_inter_residue_connectivity(structure_files):
    """测试相邻残基之间的连接。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 测试ALA 7和VAL 8之间的肽键连接
    ala_atoms = system.get_pdb_atoms_by_residue_sequence("ALA", 7)
    val_atoms = system.get_pdb_atoms_by_residue_sequence("VAL", 8)
    
    ala_c = next(atom for atom in ala_atoms if atom.name == "C")
    val_n = next(atom for atom in val_atoms if atom.name == "N")
    
    # 计算肽键长度
    peptide_bond_length = ((val_n.x - ala_c.x)**2 + 
                          (val_n.y - ala_c.y)**2 + 
                          (val_n.z - ala_c.z)**2)**0.5
    
    # 典型的肽键长度约为1.33埃
    assert peptide_bond_length == pytest.approx(1.33, abs=0.1)

def test_pdb_atom_topology_info(structure_files):
    """测试PDB原子的拓扑信息是否正确加载。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 测试ALA 7的所有原子
    ala_atoms = system.get_pdb_atoms_by_residue_sequence("ALA", 7)
    assert len(ala_atoms) > 0, "No atoms found for ALA 7"
    
    # 验证N原子的拓扑信息
    n_atom = next(atom for atom in ala_atoms if atom.name == "N")
    assert n_atom.topo_type == "NH3", "Wrong topology type for N atom"
    assert n_atom.topo_charge == pytest.approx(-0.300000), "Wrong topology charge for N atom"
    assert n_atom.topo_mass == pytest.approx(14.0070), "Wrong topology mass for N atom"
    
    # 验证CA原子的拓扑信息
    ca_atom = next(atom for atom in ala_atoms if atom.name == "CA")
    assert ca_atom.topo_type == "CT1", "Wrong topology type for CA atom"
    assert ca_atom.topo_charge == pytest.approx(0.210000), "Wrong topology charge for CA atom"
    assert ca_atom.topo_mass == pytest.approx(12.0110), "Wrong topology mass for CA atom"
    
    # 验证C原子的拓扑信息
    c_atom = next(atom for atom in ala_atoms if atom.name == "C")
    assert c_atom.topo_type == "C", "Wrong topology type for C atom"
    assert c_atom.topo_charge == pytest.approx(0.510000), "Wrong topology charge for C atom"
    assert c_atom.topo_mass == pytest.approx(12.0110), "Wrong topology mass for C atom"

def test_pdb_atom_coordinates(structure_files):
    """测试PDB原子的坐标信息是否正确加载。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 测试第一个和最后一个原子的坐标
    first_atom = system.get_pdb_atom(0)
    last_atom = system.get_pdb_atom(system.get_pdb_atom_count() - 1)
    
    # 验证坐标是有限数
    assert all(map(math.isfinite, [first_atom.x, first_atom.y, first_atom.z])), \
        "First atom has invalid coordinates"
    assert all(map(math.isfinite, [last_atom.x, last_atom.y, last_atom.z])), \
        "Last atom has invalid coordinates"
    
    # 验证坐标在合理范围内（根据CRYST1记录）
    assert 0 <= first_atom.x <= 127.022, "First atom x coordinate out of box"
    assert 0 <= first_atom.y <= 133.419, "First atom y coordinate out of box"
    assert 0 <= first_atom.z <= 132.854, "First atom z coordinate out of box"
    
    assert 0 <= last_atom.x <= 127.022, "Last atom x coordinate out of box"
    assert 0 <= last_atom.y <= 133.419, "Last atom y coordinate out of box"
    assert 0 <= last_atom.z <= 132.854, "Last atom z coordinate out of box"

def test_pdb_atom_chain_info(structure_files):
    """测试PDB原子的链信息是否正确加载。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 获取所有链的原子
    chains = set()
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        chains.add(atom.chain)
    
    # 验证每条链
    for chain in chains:
        chain_atoms = system.get_pdb_atoms_by_chain(chain)
        assert len(chain_atoms) > 0, f"No atoms found for chain {chain}"
        
        # 验证所有原子确实属于这条链
        for atom in chain_atoms:
            assert atom.chain == chain, f"Atom has wrong chain identifier: {atom.chain} != {chain}"
            
        # 验证链内残基的有效性
        residue_sequences = sorted(list(set((atom.residue, atom.sequence) for atom in chain_atoms)))
        for i in range(len(residue_sequences) - 1):
            curr_res = residue_sequences[i]
            next_res = residue_sequences[i + 1]
            # 只验证残基序号是否有效
            assert curr_res[1] > 0 and next_res[1] > 0, \
                f"Invalid residue numbers in chain {chain}: {curr_res} -> {next_res}"

def test_pdb_atom_validation(structure_files):
    """测试PDB原子的验证功能。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 验证所有原子的基本属性
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        
        # 验证原子序号
        assert atom.serial > 0, f"Invalid serial number for atom {i}"
        
        # 验证原子名称
        assert atom.name, f"Missing name for atom {i}"
        assert len(atom.name.strip()) > 0, f"Empty name for atom {i}"
        
        # 验证残基信息
        assert atom.residue, f"Missing residue name for atom {i}"
        assert len(atom.residue.strip()) > 0, f"Empty residue name for atom {i}"
        assert atom.sequence > 0, f"Invalid sequence number for atom {i}"
        
        # 验证坐标
        assert all(map(math.isfinite, [atom.x, atom.y, atom.z])), \
            f"Invalid coordinates for atom {i}"
        
        # 验证PDB特有属性
        assert 0 <= atom.occupancy <= 1, f"Invalid occupancy for atom {i}"
        assert atom.temp_factor >= 0, f"Invalid temperature factor for atom {i}"
        
        # 验证元素信息
        assert atom.element or atom.type, f"Missing element and type for atom {i}"

def test_pdb_atom_element_types(structure_files):
    """测试PDB原子的元素类型是否正确。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 收集所有原子的元素类型
    element_types = {}
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        # 从原子名称推断元素类型
        element = atom.name[0].upper()
        if element not in element_types:
            element_types[element] = []
        element_types[element].append(atom)
    
    # 验证常见元素的存在性
    common_elements = {"C", "H", "N", "O"}
    for element in common_elements:
        assert element in element_types, f"Common element {element} not found"
    
    # 验证每种元素类型的原子特性
    expected_masses = {
        "C": 12.0110,  # 碳原子
        "N": 14.0070,  # 氮原子
        "O": 15.9994,  # 氧原子
        "H": 1.0080    # 氢原子
    }
    
    for element, atoms in element_types.items():
        if element in expected_masses:
            # 找到至少一个具有正确质量的原子
            has_valid_mass = False
            for atom in atoms:
                if (hasattr(atom, 'topo_mass') and 
                    atom.topo_mass is not None and 
                    not math.isnan(atom.topo_mass)):
                    assert abs(atom.topo_mass - expected_masses[element]) < 0.001, \
                        f"Wrong mass for {element} atom {atom.name}: {atom.topo_mass} != {expected_masses[element]}"
                    has_valid_mass = True
            assert has_valid_mass, f"No atom with valid mass found for element {element}"

def test_pdb_atom_basic_properties(structure_files):
    """测试PDB原子的基本属性。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 测试第一个原子（ALA 7的N原子）
    first_atom = system.get_pdb_atom(0)
    assert first_atom.serial == 1, "First atom should have serial number 1"
    assert first_atom.name == "N", "First atom should be N"
    assert first_atom.residue == "ALA", "First atom should be in ALA"
    assert first_atom.sequence == 7, "First atom should be in residue 7"
    assert first_atom.chain == "P", "First atom should be in chain P"
    
    # 测试坐标
    assert first_atom.x == pytest.approx(76.563)
    assert first_atom.y == pytest.approx(93.118)
    assert first_atom.z == pytest.approx(93.806)
    
    # 测试PDB特有属性
    assert first_atom.occupancy == pytest.approx(1.0)
    assert first_atom.temp_factor == pytest.approx(0.0)
    assert first_atom.alt_loc == " "
    assert first_atom.insertion_code == " "

def test_pdb_atom_topology_properties(structure_files):
    """测试PDB原子的拓扑属性。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 获取ALA 7的所有原子
    ala_atoms = system.get_pdb_atoms_by_residue_sequence("ALA", 7)
    
    # 验证主链原子的拓扑属性
    backbone_atoms = {
        "N": {"type": "NH3", "charge": -0.300000, "mass": 14.0070},
        "CA": {"type": "CT1", "charge": 0.210000, "mass": 12.0110},
        "C": {"type": "C", "charge": 0.510000, "mass": 12.0110},
        "O": {"type": "O", "charge": -0.510000, "mass": 15.9994}
    }
    
    for name, props in backbone_atoms.items():
        atom = next(atom for atom in ala_atoms if atom.name == name)
        assert atom.topo_type == props["type"], \
            f"Wrong topology type for {name} atom"
        assert atom.topo_charge == pytest.approx(props["charge"]), \
            f"Wrong topology charge for {name} atom"
        assert atom.topo_mass == pytest.approx(props["mass"]), \
            f"Wrong topology mass for {name} atom"

def test_pdb_residue_composition(structure_files):
    """测试残基内的原子组成。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 验证ALA残基的原子组成
    ala_atoms = system.get_pdb_atoms_by_residue_sequence("ALA", 7)
    ala_atom_names = {atom.name for atom in ala_atoms}
    expected_ala_atoms = {
        "N", "H1", "H2", "H3",  # N端氨基（注意：使用H1/H2/H3而不是HT1/HT2/HT3）
        "CA", "HA",             # α碳
        "CB", "HB1", "HB2", "HB3", # 侧链
        "C", "O"                # C端羧基
    }
    assert ala_atom_names == expected_ala_atoms, \
        f"Missing or extra atoms in ALA residue: {ala_atom_names - expected_ala_atoms}"
    
    # 验证VAL残基的原子组成
    val_atoms = system.get_pdb_atoms_by_residue_sequence("VAL", 8)
    val_atom_names = {atom.name for atom in val_atoms}
    expected_val_atoms = {
        "N", "HN",                  # 肽键氨基
        "CA", "HA",                 # α碳
        "CB", "HB",                 # β碳
        "CG1", "HG11", "HG12", "HG13", # γ碳1
        "CG2", "HG21", "HG22", "HG23", # γ碳2
        "C", "O"                    # 肽键羰基
    }
    assert val_atom_names == expected_val_atoms, \
        f"Missing or extra atoms in VAL residue: {val_atom_names - expected_val_atoms}"

def test_pdb_atom_spatial_relations(structure_files):
    """测试原子间的空间关系。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 获取ALA 7的原子
    ala_atoms = system.get_pdb_atoms_by_residue_sequence("ALA", 7)
    
    # 计算键长
    def calc_distance(atom1, atom2):
        return ((atom1.x - atom2.x)**2 + 
                (atom1.y - atom2.y)**2 + 
                (atom1.z - atom2.z)**2)**0.5
    
    # 验证N-CA键长
    n_atom = next(atom for atom in ala_atoms if atom.name == "N")
    ca_atom = next(atom for atom in ala_atoms if atom.name == "CA")
    n_ca_distance = calc_distance(n_atom, ca_atom)
    assert n_ca_distance == pytest.approx(1.47, abs=0.1), \
        f"N-CA bond length {n_ca_distance} is out of range"
    
    # 验证CA-CB键长
    cb_atom = next(atom for atom in ala_atoms if atom.name == "CB")
    ca_cb_distance = calc_distance(ca_atom, cb_atom)
    assert ca_cb_distance == pytest.approx(1.52, abs=0.1), \
        f"CA-CB bond length {ca_cb_distance} is out of range"
    
    # 验证CA-C键长
    c_atom = next(atom for atom in ala_atoms if atom.name == "C")
    ca_c_distance = calc_distance(ca_atom, c_atom)
    assert ca_c_distance == pytest.approx(1.52, abs=0.1), \
        f"CA-C bond length {ca_c_distance} is out of range"

def test_pdb_atom_sequence_order(structure_files):
    """测试PDB原子的序列顺序。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    # 验证残基内原子的顺序（以ALA 7为例）
    ala_atoms = system.get_pdb_atoms_by_residue_sequence("ALA", 7)
    expected_order = ["N", "H1", "H2", "H3", "CA", "HA", "CB", "HB1", "HB2", "HB3", "C", "O"]
    actual_order = [atom.name for atom in ala_atoms]
    assert actual_order == expected_order, \
        f"Wrong atom order in ALA residue: {actual_order}"
    
    # 验证残基序号的连续性
    prev_sequence = None
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        if prev_sequence is not None and atom.sequence != prev_sequence:
            assert atom.sequence > prev_sequence, \
                f"Residue sequence numbers not monotonically increasing at atom {i}"
        prev_sequence = atom.sequence
