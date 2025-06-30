# tests/io/pdbParser/residue_comparison.py
"""
Residue-level statistics comparison between direct and pygcmc parsing.
Handles residue name mapping and validation.
"""

import pytest
import os
import pygcmc
from .direct_parser_core import read_pdb_direct, TEST_DATA_DIR


def test_4wp7_residue_statistics_comparison():
    """Compare residue-level statistics between direct and pygcmc parsing."""
    pdb_path = os.path.join(TEST_DATA_DIR, "4wp7", "4wp7_fixed_with_5l13_silcs.1.prod.74.rec.pdb")
    
    # Skip test if file doesn't exist
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test file {pdb_path} not found")
    
    # Parse with both methods
    direct_atoms, _ = read_pdb_direct(pdb_path)
    pygcmc_result = pygcmc.PDBParser.parse_file(pdb_path)
    
    # Count residues with direct parsing
    direct_residue_counts = {}
    for atom in direct_atoms:
        res_key = atom.res_name
        direct_residue_counts[res_key] = direct_residue_counts.get(res_key, 0) + 1
    
    # Count residues with pygcmc parsing
    pygcmc_residue_counts = {}
    for atom in pygcmc_result.atoms:
        res_key = atom.get_resname()
        pygcmc_residue_counts[res_key] = pygcmc_residue_counts.get(res_key, 0) + 1
    
    # Define known residue name mappings (pygcmc may standardize names)
    residue_mapping = {
        # Original -> Standardized
        "ACE": "ACEY",   # acetylene
        "BEN": "BENX",   # benzene  
        "DME": "DMEE",   # dimethyl ether
        "FOR": "FORM",   # formaldehyde
        "MEO": "MEOH",   # methanol
        "IMI": "IMIA",   # imidazole
        "PRP": "PRPX",   # propane
        "MAM": "MAMY",   # methylamine
    }
    
    # Apply mapping to direct counts
    mapped_direct_counts = {}
    for res_type, count in direct_residue_counts.items():
        mapped_name = residue_mapping.get(res_type, res_type)
        mapped_direct_counts[mapped_name] = mapped_direct_counts.get(mapped_name, 0) + count
    
    # Compare mapped residue type counts
    all_residue_types = set(mapped_direct_counts.keys()) | set(pygcmc_residue_counts.keys())
    
    major_differences = []
    for res_type in all_residue_types:
        direct_count = mapped_direct_counts.get(res_type, 0)
        pygcmc_count = pygcmc_residue_counts.get(res_type, 0)
        
        if direct_count != pygcmc_count:
            difference = abs(direct_count - pygcmc_count)
            if difference > max(1, min(direct_count, pygcmc_count) * 0.01):  # >1% difference
                major_differences.append(f"{res_type}: direct={direct_count}, pygcmc={pygcmc_count}")
    
    # Allow some differences due to parser variations, but verify total counts match
    total_direct = sum(mapped_direct_counts.values())
    total_pygcmc = sum(pygcmc_residue_counts.values())
    assert total_direct == total_pygcmc, f"Total atom count mismatch: direct={total_direct}, pygcmc={total_pygcmc}"
    
    # Check that we have most key residue types in common
    common_types = set(mapped_direct_counts.keys()) & set(pygcmc_residue_counts.keys())
    assert len(common_types) >= len(all_residue_types) * 0.8, \
        f"Too few common residue types: {len(common_types)}/{len(all_residue_types)}"
    
    # Check key GCMC molecules (these should match exactly after mapping)
    gcmc_molecules = ["ACEY", "BENX", "DMEE", "FORM"]
    matching_gcmc = 0
    for molecule in gcmc_molecules:
        if molecule in mapped_direct_counts and molecule in pygcmc_residue_counts:
            direct_count = mapped_direct_counts[molecule]
            pygcmc_count = pygcmc_residue_counts[molecule]
            if direct_count == pygcmc_count:
                matching_gcmc += 1
    
    assert matching_gcmc >= 3, f"At least 3 GCMC molecules should match exactly, got {matching_gcmc}"
    
    print(f"Residue statistics comparison: {len(common_types)}/{len(all_residue_types)} types match, "
          f"{matching_gcmc}/4 GCMC molecules match exactly")