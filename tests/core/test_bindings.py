# tests/core/test_bindings.py

import pytest
import pygcmc
import math

@pytest.fixture
def water_system():
    """创建一个包含单个水分子的系统"""
    # 创建水分子的三个原子
    h1 = pygcmc.PDBAtom()
    h1.serial = 1
    h1.name = "H1"
    h1.residue = "HOH"
    h1.sequence = 1
    h1.chain = "A"
    h1.x = 0.0
    h1.y = 0.0
    h1.z = 0.0
    h1.element = "H"
    h1.type = "H"
    h1.topo_type = "HT"
    h1.topo_charge = 0.417
    h1.topo_mass = 1.008
    h1.forcefield_epsilon = 0.0157
    h1.forcefield_rmin = 0.6000
    h1.occupancy = 1.0
    h1.temp_factor = 0.0
    h1.alt_loc = " "
    h1.insertion_code = " "

    o1 = pygcmc.PDBAtom()
    o1.serial = 2
    o1.name = "O"
    o1.residue = "HOH"
    o1.sequence = 1
    o1.chain = "A"
    o1.x = 1.0
    o1.y = 0.0
    o1.z = 0.0
    o1.element = "O"
    o1.type = "O"
    o1.topo_type = "OT"
    o1.topo_charge = -0.834
    o1.topo_mass = 15.999
    o1.forcefield_epsilon = 0.1521
    o1.forcefield_rmin = 1.7682
    o1.occupancy = 1.0
    o1.temp_factor = 0.0
    o1.alt_loc = " "
    o1.insertion_code = " "

    h2 = pygcmc.PDBAtom()
    h2.serial = 3
    h2.name = "H2"
    h2.residue = "HOH"
    h2.sequence = 1
    h2.chain = "A"
    h2.x = 1.0
    h2.y = 1.0
    h2.z = 0.0
    h2.element = "H"
    h2.type = "H"
    h2.topo_type = "HT"
    h2.topo_charge = 0.417
    h2.topo_mass = 1.008
    h2.forcefield_epsilon = 0.0157
    h2.forcefield_rmin = 0.6000
    h2.occupancy = 1.0
    h2.temp_factor = 0.0
    h2.alt_loc = " "
    h2.insertion_code = " "
    
    # 创建水分子残基
    water = pygcmc.IOResidue()
    water.name = "HOH"
    water.sequence_number = 1
    water.chain_id = 'A'
    water.atoms = [h1, o1, h2]
    
    return water, [h1, o1, h2]

def test_pdb_atom_properties():
    """测试PDBAtom的属性设置和获取"""
    atom = pygcmc.PDBAtom()
    # 基本属性
    atom.serial = 1
    atom.name = "H1"
    atom.residue = "HOH"
    atom.sequence = 1
    atom.chain = "A"
    atom.x = 1.0
    atom.y = 2.0
    atom.z = 3.0
    atom.element = "H"
    atom.type = "H"
    # 拓扑属性
    atom.topo_type = "HT"
    atom.topo_charge = 0.417
    atom.topo_mass = 1.008
    # 力场属性
    atom.forcefield_epsilon = 0.0157
    atom.forcefield_rmin = 0.6000
    # PDB特有属性
    atom.occupancy = 1.0
    atom.temp_factor = 0.0
    atom.alt_loc = " "
    atom.insertion_code = " "
    
    # 验证所有属性
    assert atom.serial == 1
    assert atom.name == "H1"
    assert atom.residue == "HOH"
    assert atom.sequence == 1
    assert atom.chain == "A"
    assert atom.x == 1.0
    assert atom.y == 2.0
    assert atom.z == 3.0
    assert atom.element == "H"
    assert atom.type == "H"
    assert atom.topo_type == "HT"
    assert abs(atom.topo_charge - 0.417) < 1e-6
    assert abs(atom.topo_mass - 1.008) < 1e-6
    assert abs(atom.forcefield_epsilon - 0.0157) < 1e-6
    assert abs(atom.forcefield_rmin - 0.6000) < 1e-6
    assert atom.occupancy == 1.0
    assert atom.temp_factor == 0.0
    assert atom.alt_loc == " "
    assert atom.insertion_code == " "
    assert atom.is_valid()
    assert atom.has_topology_info()
    assert atom.has_forcefield_info()

def test_residue_properties(water_system):
    """测试IOResidue的属性和方法"""
    residue, atoms = water_system
    # 基本属性测试
    assert residue.name == "HOH"
    assert residue.sequence_number == 1
    assert residue.chain_id == 'A'
    assert len(residue.atoms) == 3
    
    # 测试center_of_mass方法
    com = residue.center_of_mass()
    assert len(com) == 3
    assert isinstance(com[0], float)
    assert isinstance(com[1], float)
    assert isinstance(com[2], float)
    
    # 测试atom_count方法
    assert residue.atom_count() == 3
    
    # 测试原子列表访问
    assert all(isinstance(atom, pygcmc.PDBAtom) for atom in residue.atoms)
    assert residue.atoms[0].name == "H1"
    assert residue.atoms[1].name == "O"
    assert residue.atoms[2].name == "H2"

def test_forcefield_pair():
    """测试ForceFieldPair的属性和方法"""
    # 测试默认构造函数
    ff_pair = pygcmc.ForceFieldPair()
    assert hasattr(ff_pair, 'rmin')
    assert hasattr(ff_pair, 'epsilon')
    
    # 测试带参数的构造函数
    ff_pair2 = pygcmc.ForceFieldPair(1.0, 2.0)
    assert abs(ff_pair2.rmin - 1.0) < 1e-6
    assert abs(ff_pair2.epsilon - 2.0) < 1e-6
    
    # 测试属性修改
    ff_pair2.rmin = 1.5
    ff_pair2.epsilon = 2.5
    assert abs(ff_pair2.rmin - 1.5) < 1e-6
    assert abs(ff_pair2.epsilon - 2.5) < 1e-6

def test_project_atom():
    """测试ProjectAtom的属性和方法"""
    # 从PDBAtom创建ProjectAtom
    pdb_atom = pygcmc.PDBAtom()
    pdb_atom.serial = 1
    pdb_atom.name = "H1"
    pdb_atom.residue = "HOH"
    pdb_atom.sequence = 1
    pdb_atom.chain = "A"
    pdb_atom.x = 1.0
    pdb_atom.y = 2.0
    pdb_atom.z = 3.0
    
    project_atom = pygcmc.ProjectAtom(pdb_atom)
    
    # 验证属性复制
    assert project_atom.serial == 1
    assert project_atom.name == "H1"
    assert project_atom.residue == "HOH"
    assert project_atom.sequence == 1
    assert project_atom.chain == "A"
    assert project_atom.x == 1.0
    assert project_atom.y == 2.0
    assert project_atom.z == 3.0

def test_pdb_parser():
    """测试PDBParser的功能"""
    parser = pygcmc.PDBParser()
    # 注意：这里需要一个实际的PDB文件来测试parse方法
    # 目前只测试实例化

def test_top_parser():
    """测试TopParser的功能"""
    parser = pygcmc.TopParser()
    # 测试基本实例化
    assert hasattr(parser, 'parse')
    assert hasattr(parser, 'parse_with_includes')
    assert hasattr(parser, 'update_pdb_atoms')

def test_ff_parser():
    """测试FFParser的功能"""
    parser = pygcmc.FFParser()
    # 测试基本属性和方法
    assert hasattr(parser, 'parse')
    assert hasattr(parser, 'get_nonbonded_params')
    assert hasattr(parser, 'get_nbfix_params')
    assert hasattr(parser, 'get_cutnb')
    assert hasattr(parser, 'get_ctofnb')
    assert hasattr(parser, 'get_ctonnb')
    assert hasattr(parser, 'get_eps')
    assert hasattr(parser, 'get_e14fac')
    assert hasattr(parser, 'get_wmin')