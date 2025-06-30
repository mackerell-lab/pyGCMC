# tests/io/topParser/enhanced_4wp7_advanced.py
"""
Enhanced 4wp7 TOP Parser Tests - Advanced Features

This file contains advanced feature tests for the 4wp7 topology file.

Tests include:
- Preprocessor directive handling
- Error handling and edge cases
"""

import pytest
import os
import tempfile
from pygcmc.io import TOPParser
from pygcmc.model import Topology

# Get the directory containing test data files
TEST_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "data")


def test_4wp7_preprocessor_directives():
    """Test preprocessor directive handling (#ifdef, #include, #define)."""
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    # Read the raw file to check preprocessor directives
    with open(top_path, 'r') as f:
        content = f.read()
    
    # Verify key preprocessor features are present
    assert "#include" in content, "Should contain #include directives"
    
    # Count include directives
    include_count = content.count("#include")
    assert include_count >= 10, f"Expected at least 10 includes (forcefield + molecules + water), got {include_count}"
    
    # Verify force field include is first
    lines = content.split('\n')
    first_include_line = None
    for line in lines:
        if line.strip().startswith("#include"):
            first_include_line = line.strip()
            break
    
    assert first_include_line is not None, "Should have at least one #include directive"
    assert "forcefield.itp" in first_include_line, f"First include should be forcefield, got: {first_include_line}"
    
    # Parse with the topology parser to verify preprocessor handling
    parser = TOPParser()
    topology = Topology()
    
    success = parser.parse_to_topology(top_path, topology)
    assert success, "Failed to parse 4wp7 TOP file with preprocessor directives"
    
    # Verify the parser correctly processed includes by checking for water molecules
    # Water should only be present if tip3p.itp was correctly processed with _FF_CHARMM defined
    water_found = False
    for i in range(topology.get_num_residues()):
        residue = topology.get_residue(i)
        if residue.name == "SOL":
            water_found = True
            break
    
    assert water_found, "Water molecules not found - preprocessor may not be handling _FF_CHARMM correctly"
    
    # Verify all expected includes were processed by checking molecule presence
    expected_molecules = {"BENX", "PRPX", "DMEE", "MEOH", "FORM", "IMIA", "ACEY", "MAMY", "SOL"}
    found_molecules = set()
    
    for i in range(topology.get_num_residues()):
        residue = topology.get_residue(i)
        if residue.name in expected_molecules:
            found_molecules.add(residue.name)
    
    missing_molecules = expected_molecules - found_molecules
    assert len(missing_molecules) == 0, \
        f"Missing molecules suggest include processing failed: {missing_molecules}"
    
    print(f"✓ Preprocessor directives validated: {include_count} includes processed, all molecules found")


def test_4wp7_error_handling():
    """Test error handling for corrupted or incomplete 4wp7-style files."""
    top_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.gc.74.top")
    
    # Skip test if file doesn't exist
    if not os.path.exists(top_path):
        pytest.skip(f"Test file {top_path} not found")
    
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Test 1: Missing include file
        corrupted_file1 = os.path.join(tmp_dir, "missing_include.top")
        with open(corrupted_file1, 'w') as f:
            f.write('#include "./nonexistent.itp"\n')
            f.write('[ moleculetype ]\n')
            f.write('Test 1\n')
            f.write('[ atoms ]\n')
            f.write('1 CT1 1 TEST C1 1 0.0 12.011\n')
            f.write('[ system ]\n')
            f.write('Test\n')
            f.write('[ molecules ]\n')
            f.write('Test 1\n')
        
        parser = TOPParser()
        topology = Topology()
        
        # Should handle missing include gracefully (may succeed or fail, but shouldn't crash)
        try:
            result = parser.parse_to_topology(corrupted_file1, topology)
            # Either succeeds (graceful handling) or fails (proper error)
            assert isinstance(result, bool), "Parser should return boolean result"
        except Exception as e:
            # If it throws exception, should be a reasonable error message
            assert len(str(e)) > 0, "Exception should have meaningful message"
        
        # Test 2: Malformed [molecules] section
        corrupted_file2 = os.path.join(tmp_dir, "bad_molecules.top")
        with open(corrupted_file2, 'w') as f:
            f.write('[ moleculetype ]\n')
            f.write('Test 1\n')
            f.write('[ atoms ]\n')
            f.write('1 CT1 1 TEST C1 1 0.0 12.011\n')
            f.write('[ system ]\n')
            f.write('Test\n')
            f.write('[ molecules ]\n')
            f.write('Test 999999\n')  # Impossibly high count
        
        topology2 = Topology()
        
        # Should handle unrealistic molecule counts gracefully
        try:
            result = parser.parse_to_topology(corrupted_file2, topology2)
            if result:
                # If parsing succeeds, should not create impossibly many atoms
                assert topology2.get_num_atoms() <= 1000000, "Should not create excessive atoms from bad molecule count"
        except Exception as e:
            # Acceptable to fail with reasonable error
            assert "molecule" in str(e).lower() or "count" in str(e).lower(), \
                f"Error message should mention molecule/count issue: {e}"
        
        # Test 3: Truncated file (incomplete atoms section)
        truncated_file = os.path.join(tmp_dir, "truncated.top")
        with open(top_path, 'r') as original:
            lines = original.readlines()
        
        # Take first 1000 lines (cuts off in middle of atoms section)
        with open(truncated_file, 'w') as f:
            f.writelines(lines[:1000])
        
        topology3 = Topology()
        result = parser.parse_to_topology(truncated_file, topology3)
        
        # Truncated file might succeed (partial parsing) or fail, but shouldn't crash
        assert isinstance(result, bool), "Parser should handle truncated file gracefully"
        
        if result:
            # If it succeeds, should have parsed some atoms but not the full system
            assert 0 < topology3.get_num_atoms() < 50000, \
                f"Truncated parsing should give partial result, got {topology3.get_num_atoms()} atoms"
    
    print("✓ Error handling validated: missing includes, bad molecule counts, truncated files")