"""
Test ActivePool memory management
"""

import pytest
from .fixtures import (
    MOVEMENT_AVAILABLE,
    create_active_pool,
    create_water_molecule
)


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_pool_creation():
    """Test pool creation with capacity"""
    active_pool = create_active_pool()
    assert active_pool.getMaxAtoms() == 1000
    assert active_pool.getMaxResidues() == 100
    counts = active_pool.getActiveCounts()
    assert counts == (0, 0)  # Initially empty


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_insert_molecule():
    """Test molecule insertion into pool"""
    active_pool = create_active_pool()
    atoms = create_water_molecule()
    res_idx = active_pool.insertMolecule(atoms, resType=0)
    assert res_idx >= 0
    assert active_pool.isResidueActive(res_idx)
    
    counts = active_pool.getActiveCounts()
    assert counts[0] == len(atoms)  # Active atoms
    assert counts[1] == 1  # Active residues


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_delete_residue():
    """Test residue deletion from pool"""
    active_pool = create_active_pool()
    atoms = create_water_molecule()
    res_idx = active_pool.insertMolecule(atoms, resType=0)
    
    success = active_pool.deleteResidue(res_idx)
    assert success
    assert not active_pool.isResidueActive(res_idx)
    
    counts = active_pool.getActiveCounts()
    assert counts == (0, 0)  # Should be empty after deletion


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_fragmentation():
    """Test fragmentation calculation and compaction"""
    active_pool = create_active_pool()
    
    # Insert and delete to create fragmentation
    indices = []
    for _ in range(10):
        atoms = create_water_molecule()
        idx = active_pool.insertMolecule(atoms, resType=0)
        indices.append(idx)
    
    # Delete every other molecule
    for i in range(0, 10, 2):
        active_pool.deleteResidue(indices[i])
    
    fragmentation = active_pool.getFragmentation()
    assert 0.0 <= fragmentation <= 1.0
    
    # Compact the pool
    compacted = active_pool.compact(force=True)
    assert compacted >= 0  # May or may not compact depending on fragmentation


@pytest.mark.skipif(not MOVEMENT_AVAILABLE, reason="Movement module not compiled")
def test_capacity_check():
    """Test capacity checking"""
    active_pool = create_active_pool()
    can_insert = active_pool.canInsert(3)  # 3 atoms for water
    assert can_insert  # Should have space initially
    
    # Fill the pool (but not completely to avoid infinite loop)
    for _ in range(100):  # Insert up to 100 molecules
        if active_pool.canInsert(3):
            atoms = create_water_molecule()
            active_pool.insertMolecule(atoms, resType=0)
        else:
            break
    
    # Check that we inserted some molecules
    counts = active_pool.getActiveCounts()
    assert counts[1] > 0  # Should have some residues