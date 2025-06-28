# tests/system/molecularSystem/incompatible_tests.py

import pytest
from .helpers import *

def test_combine_incompatible_files():
    """Test error handling when combining incompatible structure and topology files."""
    # Load test.pdb which contains protein, BENX, PRPX, and SOL
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Load test_proa.psf which only contains protein chain
    psf_path = os.path.join(TEST_DATA_DIR, "test_proa.psf")
    topology = pygcmc.io.PSFParser.parse_file(psf_path)
    
    # Attempt to combine should raise an error
    mol_system = pygcmc.MolecularSystem()
    with pytest.raises(RuntimeError) as excinfo:
        molecular = mol_system.combine(structure, topology)
    
    # Check that the error message is descriptive
    error_msg = str(excinfo.value)
    assert "Inconsistent total number of" in error_msg
    assert "Structure has" in error_msg
    assert "but Topology has" in error_msg

def test_combine_multiple_topologies():
    """Test combining structure with multiple topology files."""
    # Load structure file
    pdb_path = os.path.join(TEST_DATA_DIR, "test.pdb")
    structure = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Load topology files
    proa_path = os.path.join(TEST_DATA_DIR, "test_proa.psf")
    benx_path = os.path.join(TEST_DATA_DIR, "mols/benx.psf")
    prpx_path = os.path.join(TEST_DATA_DIR, "mols/prpx.itp")
    sol_path = os.path.join(TEST_DATA_DIR, "mols/sol.itp")
    
    topologies = [
        pygcmc.io.PSFParser.parse_file(proa_path),
        pygcmc.io.PSFParser.parse_file(benx_path),
        pygcmc.io.TOPParser.parse_file(prpx_path),
        pygcmc.io.TOPParser.parse_file(sol_path)
    ]
    
    # Combine structure with multiple topologies
    mol_system = pygcmc.MolecularSystem()
    molecular = mol_system.combine_multiple(structure, topologies)
    
    # Verify the combined result
    assert len(molecular.atoms) == len(structure.atoms)
    assert len(molecular.residues) == len(structure.residues)
    
    # Check specific residues
    residues = molecular.residues
    assert residues[0].get_resname() == "ALA"  # From test_proa.psf
    assert residues[1].get_resname() == "VAL"  # From test_proa.psf
    assert residues[2].get_resname() == "PRO"  # From test_proa.psf
    
    # Find BENX residue
    benx_found = False
    for res in residues:
        if res.get_resname() == "BENX":
            benx_found = True
            break
    assert benx_found
    
    # Find PRPX residues
    prpx_count = sum(1 for res in residues if res.get_resname() == "PRPX")
    assert prpx_count == 2
    
    # Find SOL residues
    sol_count = sum(1 for res in residues if res.get_resname() == "SOL")
    assert sol_count == 10

