# tests/core/test_system_multi_psf.py

import os
import pytest
import math
from pygcmc import System

@pytest.fixture
def test_data_dir():
    """获取测试数据目录的路径。"""
    current_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(current_dir, "..", "data")

def test_load_structure_psf_solvent(test_data_dir):
    """测试使用sol.psf加载溶剂分子的情况。"""
    # 准备文件路径
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    psf_file = os.path.join(test_data_dir, "mols", "sol.psf")
    
    # 加载系统
    system = System(pdb=pdb_file, psf=psf_file)
    
    # 验证系统成功加载
    assert system.has_pdb_atoms(), "System should have PDB atoms"
    
    # 收集所有溶剂分子
    solvent_atoms = []
    solvent_residues = set()
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        if atom.residue in {"SOL", "WAT", "HOH", "TIP3"}:
            solvent_atoms.append(atom)
            solvent_residues.add((atom.residue, atom.sequence))
    
    # 验证溶剂分子的存在性
    assert len(solvent_atoms) > 0, "System should have solvent atoms"
    assert len(solvent_residues) > 0, "System should have solvent residues"
    
    # 验证每个溶剂残基的完整性
    for residue_name, sequence in sorted(solvent_residues):
        atoms = system.get_pdb_atoms_by_residue_sequence(residue_name, sequence)
        assert len(atoms) == 3, f"Solvent residue {residue_name} {sequence} should have 3 atoms"
        
        # 验证原子名称（应该是 OW、HW1、HW2）
        atom_names = {atom.name for atom in atoms}
        assert "OW" in atom_names, f"Solvent residue {residue_name} {sequence} should have OW atom"
        assert "HW1" in atom_names, f"Solvent residue {residue_name} {sequence} should have HW1 atom"
        assert "HW2" in atom_names, f"Solvent residue {residue_name} {sequence} should have HW2 atom"
        
        # 验证拓扑信息
        for atom in atoms:
            assert hasattr(atom, 'topo_type'), f"Atom {atom.name} should have topology type"
            assert hasattr(atom, 'topo_charge'), f"Atom {atom.name} should have topology charge"
            assert hasattr(atom, 'topo_mass'), f"Atom {atom.name} should have topology mass"
            
            # 验证具体的拓扑信息
            if atom.name == "OW":
                assert atom.topo_type == "OT", f"Wrong topology type for OW atom"
                assert atom.topo_charge == pytest.approx(-0.834), \
                    f"Wrong topology charge for OW atom"
                assert atom.topo_mass == pytest.approx(15.9994), \
                    f"Wrong topology mass for OW atom"
            elif atom.name in {"HW1", "HW2"}:
                assert atom.topo_type == "HT", f"Wrong topology type for {atom.name} atom"
                assert atom.topo_charge == pytest.approx(0.417), \
                    f"Wrong topology charge for {atom.name} atom"
                assert atom.topo_mass == pytest.approx(1.0080), \
                    f"Wrong topology mass for {atom.name} atom"
            
        # 验证水分子的几何构型
        ow = next(atom for atom in atoms if atom.name == "OW")
        hw1 = next(atom for atom in atoms if atom.name == "HW1")
        hw2 = next(atom for atom in atoms if atom.name == "HW2")
        
        # 计算O-H键长
        oh1_length = math.sqrt(
            (ow.x - hw1.x)**2 + (ow.y - hw1.y)**2 + (ow.z - hw1.z)**2
        )
        oh2_length = math.sqrt(
            (ow.x - hw2.x)**2 + (ow.y - hw2.y)**2 + (ow.z - hw2.z)**2
        )
        
        # 验证O-H键长（应该在0.95-1.0埃范围内）
        assert 0.95 <= oh1_length <= 1.0, \
            f"Unreasonable O-H1 bond length ({oh1_length}Å) in residue {residue_name} {sequence}"
        assert 0.95 <= oh2_length <= 1.0, \
            f"Unreasonable O-H2 bond length ({oh2_length}Å) in residue {residue_name} {sequence}"
        
        # 计算H-O-H键角
        vec1 = [hw1.x - ow.x, hw1.y - ow.y, hw1.z - ow.z]
        vec2 = [hw2.x - ow.x, hw2.y - ow.y, hw2.z - ow.z]
        dot_product = sum(a * b for a, b in zip(vec1, vec2))
        len1 = math.sqrt(sum(x * x for x in vec1))
        len2 = math.sqrt(sum(x * x for x in vec2))
        angle = math.acos(dot_product / (len1 * len2)) * 180 / math.pi
        
        # 验证H-O-H键角（应该在104-106度范围内）
        assert 104.0 <= angle <= 106.0, \
            f"Unreasonable H-O-H angle ({angle}°) in residue {residue_name} {sequence}"

def test_load_structure_psf_individual(test_data_dir):
    """测试单独加载每个PSF文件的情况。"""
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    
    # 测试蛋白质PSF (test_proa.psf)
    system = System(pdb=pdb_file, psf=os.path.join(test_data_dir, "test_proa.psf"))
    assert system.has_pdb_atoms(), "System should have PDB atoms"
    atoms = system.get_pdb_atoms_by_residue("ALA")
    assert len(atoms) > 0, "No ALA atoms found"
    
    # 打印所有 ALA 残基的原子信息
    print("\nALA residue atoms:")
    for atom in atoms:
        print(f"Atom {atom.name} in ALA {atom.sequence}:")
        print(f"  - topo_type: {atom.topo_type}")
        print(f"  - topo_charge: {atom.topo_charge}")
        print(f"  - topo_mass: {atom.topo_mass}")
        if atom.name == "N":
            # Check for correct N atom type based on residue position
            if atom.sequence == 7:  # N-terminal ALA
                assert atom.topo_type == "NH3", f"Wrong topology type for N atom in N-terminal ALA"
                assert atom.topo_charge == pytest.approx(-0.30), f"Wrong topology charge for N atom in N-terminal ALA"
            else:  # Non-terminal ALA
                assert atom.topo_type == "NH1", f"Wrong topology type for N atom in non-terminal ALA"
                assert atom.topo_charge == pytest.approx(-0.47), f"Wrong topology charge for N atom in non-terminal ALA"
        elif atom.sequence == 7:  # N-terminal ALA hydrogens
            if atom.name in ["HT1", "HT2", "HT3"]:
                assert atom.topo_type == "HC", f"Wrong topology type for {atom.name} in N-terminal ALA"
                assert atom.topo_charge == pytest.approx(0.33), f"Wrong topology charge for {atom.name} in N-terminal ALA"
        else:  # Non-terminal ALA hydrogens
            if atom.name == "HN":
                assert atom.topo_type == "H", f"Wrong topology type for HN in non-terminal ALA"
                assert atom.topo_charge == pytest.approx(0.31), f"Wrong topology charge for HN in non-terminal ALA"
    
    # 测试苯分子PSF (benx.psf)
    system = System(pdb=pdb_file, psf=os.path.join(test_data_dir, "mols", "benx.psf"))
    assert system.has_pdb_atoms(), "System should have PDB atoms"
    benx_atoms = system.get_pdb_atoms_by_residue("BENX")
    assert len(benx_atoms) > 0, "No BENX atoms found"
    
    # 打印所有 BENX 残基的原子信息
    print("\nBENX residue atoms:")
    for atom in benx_atoms:
        print(f"Atom {atom.name} in BENX {atom.sequence}:")
        print(f"  - topo_type: {atom.topo_type}")
        print(f"  - topo_charge: {atom.topo_charge}")
        print(f"  - topo_mass: {atom.topo_mass}")
        if atom.name == "CG":
            assert atom.topo_type == "CG2R61", f"Wrong topology type for CG atom in BENX"
            assert atom.topo_charge == pytest.approx(-0.115), f"Wrong topology charge for CG atom in BENX"
    
    # 测试丙烷分子PSF (prpx.psf)
    system = System(pdb=pdb_file, psf=os.path.join(test_data_dir, "mols", "prpx.psf"))
    assert system.has_pdb_atoms(), "System should have PDB atoms"
    prpx_atoms = system.get_pdb_atoms_by_residue("PRPX")
    assert len(prpx_atoms) > 0, "No PRPX atoms found"
    
    # 打印所有 PRPX 残基的原子信息
    print("\nPRPX residue atoms:")
    for atom in prpx_atoms:
        print(f"Atom {atom.name} in PRPX {atom.sequence}:")
        print(f"  - topo_type: {atom.topo_type}")
        print(f"  - topo_charge: {atom.topo_charge}")
        print(f"  - topo_mass: {atom.topo_mass}")
        if atom.name == "C1":
            assert atom.topo_type == "CG331", f"Wrong topology type for C1 atom in PRPX"
            assert atom.topo_charge == pytest.approx(-0.27), f"Wrong topology charge for C1 atom in PRPX"
    
    # 测试水分子PSF (sol.psf)
    system = System(pdb=pdb_file, psf=os.path.join(test_data_dir, "mols", "sol.psf"))
    assert system.has_pdb_atoms(), "System should have PDB atoms"
    solvent_atoms = []
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        if atom.residue in {"SOL", "WAT", "HOH", "TIP3"}:
            solvent_atoms.append(atom)
    assert len(solvent_atoms) > 0, "No solvent atoms found"
    
    # 打印所有溶剂分子的原子信息
    print("\nSolvent atoms:")
    for atom in solvent_atoms:
        print(f"Atom {atom.name} in {atom.residue} {atom.sequence}:")
        print(f"  - topo_type: {atom.topo_type}")
        print(f"  - topo_charge: {atom.topo_charge}")
        print(f"  - topo_mass: {atom.topo_mass}")
        if atom.name == "OW":
            assert atom.topo_type == "OT", f"Wrong topology type for OW atom in SOL"
            assert atom.topo_charge == pytest.approx(-0.834), f"Wrong topology charge for OW atom in SOL"

def test_load_structure_psf_multi(test_data_dir):
    """测试顺序加载多个PSF文件的情况。"""
    # 准备文件路径
    pdb_file = os.path.join(test_data_dir, "test.pdb")
    psf_files = [
        os.path.join(test_data_dir, "test_proa.psf"),
        os.path.join(test_data_dir, "mols", "benx.psf"),
        os.path.join(test_data_dir, "mols", "prpx.psf"),
        os.path.join(test_data_dir, "mols", "sol.psf")
    ]
    
    # 加载系统
    system = System(pdb=pdb_file, psf=psf_files)
    
    # 验证系统成功加载
    assert system.has_pdb_atoms(), "System should have PDB atoms"
    
    # 验证蛋白质残基 (test_proa.psf)
    protein_residues = {"ALA", "VAL", "PRO", "ASN", "GLN"}
    for residue_name in protein_residues:
        atoms = system.get_pdb_atoms_by_residue(residue_name)
        assert len(atoms) > 0, f"No atoms found for protein residue {residue_name}"
        
        # 打印所有蛋白质残基的原子信息
        print(f"\n{residue_name} residue atoms:")
        for atom in atoms:
            print(f"Atom {atom.name} in {residue_name} {atom.sequence}:")
            print(f"  - topo_type: {atom.topo_type}")
            print(f"  - topo_charge: {atom.topo_charge}")
            print(f"  - topo_mass: {atom.topo_mass}")
            if residue_name == "ALA" and atom.name == "N":
                # Check for correct N atom type based on residue position
                if atom.sequence == 7:  # N-terminal ALA
                    assert atom.topo_type == "NH3", f"Wrong topology type for N atom in N-terminal ALA"
                    assert atom.topo_charge == pytest.approx(-0.30), f"Wrong topology charge for N atom in N-terminal ALA"
                else:  # Non-terminal ALA
                    assert atom.topo_type == "NH1", f"Wrong topology type for N atom in non-terminal ALA"
                    assert atom.topo_charge == pytest.approx(-0.47), f"Wrong topology charge for N atom in non-terminal ALA"
    
    # 验证苯分子 (benx.psf)
    benx_atoms = system.get_pdb_atoms_by_residue("BENX")
    assert len(benx_atoms) > 0, "No BENX atoms found"
    benx_atom_names = {atom.name for atom in benx_atoms}
    expected_benx_atoms = {"CG", "CD1", "CD2", "CE1", "CE2", "CZ", "HG", "HD1", "HD2", "HE1", "HE2", "HZ", "LPA"}
    assert benx_atom_names == expected_benx_atoms, f"Missing or extra atoms in BENX residue"
    
    # 打印所有 BENX 残基的原子信息
    print("\nBENX residue atoms:")
    for atom in benx_atoms:
        print(f"Atom {atom.name} in BENX {atom.sequence}:")
        print(f"  - topo_type: {atom.topo_type}")
        print(f"  - topo_charge: {atom.topo_charge}")
        print(f"  - topo_mass: {atom.topo_mass}")
        if atom.name == "CG":
            assert atom.topo_type == "CG2R61", f"Wrong topology type for CG atom in BENX"
            assert atom.topo_charge == pytest.approx(-0.115), f"Wrong topology charge for CG atom in BENX"
    
    # 验证丙烷分子 (prpx.psf)
    prpx_atoms = system.get_pdb_atoms_by_residue("PRPX")
    assert len(prpx_atoms) > 0, "No PRPX atoms found"
    prpx_atom_names = {atom.name for atom in prpx_atoms}
    expected_prpx_atoms = {"C1", "C2", "C3", "H11", "H12", "H13", "H21", "H22", "H31", "H32", "H33", "LPA"}
    assert prpx_atom_names == expected_prpx_atoms, f"Missing or extra atoms in PRPX residue"
    
    # 打印所有 PRPX 残基的原子信息
    print("\nPRPX residue atoms:")
    for atom in prpx_atoms:
        print(f"Atom {atom.name} in PRPX {atom.sequence}:")
        print(f"  - topo_type: {atom.topo_type}")
        print(f"  - topo_charge: {atom.topo_charge}")
        print(f"  - topo_mass: {atom.topo_mass}")
        if atom.name == "C1":
            assert atom.topo_type == "CG331", f"Wrong topology type for C1 atom in PRPX"
            assert atom.topo_charge == pytest.approx(-0.27), f"Wrong topology charge for C1 atom in PRPX"
    
    # 验证溶剂分子 (sol.psf)
    solvent_residues = set()
    for i in range(system.get_pdb_atom_count()):
        atom = system.get_pdb_atom(i)
        if atom.residue in {"SOL", "WAT", "HOH", "TIP3"}:
            solvent_residues.add((atom.residue, atom.sequence))
    
    assert len(solvent_residues) > 0, "No solvent residues found"
    
    # 验证每个溶剂残基的完整性
    for residue_name, sequence in sorted(solvent_residues):
        atoms = system.get_pdb_atoms_by_residue_sequence(residue_name, sequence)
        assert len(atoms) == 3, f"Solvent residue {residue_name} {sequence} should have 3 atoms"
        atom_names = {atom.name for atom in atoms}
        assert atom_names == {"OW", "HW1", "HW2"}, f"Wrong atom names in solvent residue {residue_name} {sequence}"
        
        # 打印所有溶剂分子的原子信息
        print(f"\nSolvent residue {residue_name} {sequence} atoms:")
        for atom in atoms:
            print(f"Atom {atom.name} in {residue_name} {sequence}:")
            print(f"  - topo_type: {atom.topo_type}")
            print(f"  - topo_charge: {atom.topo_charge}")
            print(f"  - topo_mass: {atom.topo_mass}")
            if atom.name == "OW":
                assert atom.topo_type == "OT", f"Wrong topology type for OW atom in {residue_name}"
                assert atom.topo_charge == pytest.approx(-0.834), f"Wrong topology charge for OW atom in {residue_name}"