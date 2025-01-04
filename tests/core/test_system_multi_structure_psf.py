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

def test_load_structure_psf_multi(structure_files):
    """测试使用多残基 PSF 文件加载结构。"""
    system = System()
    system.load_structure_psf_multi(structure_files['pdb'], structure_files['psf'])
    
    expected_residue_count = get_pdb_residue_count(structure_files['pdb'])
    verify_system_content(system, expected_residue_count=expected_residue_count)

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

def test_load_structure_psf_auto_multi(structure_files):
    """测试自动加载多残基 PSF 文件。"""
    system = System()
    system.load_structure_psf_auto(structure_files['pdb'], structure_files['psf'])
    
    expected_residue_count = get_pdb_residue_count(structure_files['pdb'])
    verify_system_content(system, expected_residue_count=expected_residue_count)

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
