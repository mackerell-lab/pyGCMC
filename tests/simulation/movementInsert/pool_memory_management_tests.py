# pool_memory_management_tests.py
"""
Tests for active pool memory management with insert/delete/compact operations.
Validates CPU-side implementation preparing for GPU migration.
"""

import pytest
import random
import math
import pygcmc
from .active_pool import ActivePool, BatchedActivePool, total_energy_active, calculate_active_energy_components
from .basic_insertion_helpers import create_empty_system, create_water_molecule, calculate_system_energy


def test_active_pool_basic_operations():
    """Test basic insert, delete, and compact operations."""
    # Create forcefield
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.ljSigma = [0.315, 0.0]  # O, H
    ff.ljEps = [0.636, 0.0]
    
    # Initialize pool
    pool = ActivePool(box=[5.0, 5.0, 5.0], cutoff=1.2, forcefield=ff)
    
    # Insert 3 water molecules
    res0 = pool.insert_molecule(create_water_molecule(1.0, 1.0, 1.0))
    print(f"After res0: atoms={len(pool.state.atoms)}, metadata={pool.residue_metadata}")
    res1 = pool.insert_molecule(create_water_molecule(2.0, 2.0, 2.0))
    print(f"After res1: atoms={len(pool.state.atoms)}, metadata={pool.residue_metadata}")
    res2 = pool.insert_molecule(create_water_molecule(3.0, 3.0, 3.0))
    print(f"After res2: atoms={len(pool.state.atoms)}, metadata={pool.residue_metadata}")
    
    assert pool.state.activeResidueCount == 3
    assert pool.state.activeAtomCount == 9  # 3 atoms per water
    assert pool.stats['total_inserts'] == 3
    
    # Delete middle molecule (ghost it)
    success = pool.delete_residue(res1)
    assert success
    assert pool.state.activeResidueCount == 2
    assert pool.state.activeAtomCount == 6
    assert pool.stats['total_deletes'] == 1
    
    # Check fragmentation
    frag = pool.get_fragmentation()
    assert abs(frag - 1.0/3.0) < 0.01  # 1/3 inactive
    
    # Debug info before compact
    print(f"Before compact: {len(pool.residue_metadata)} metadata entries")
    for i, m in enumerate(pool.residue_metadata):
        print(f"  {i}: active={m.active}, start={m.atom_start}, count={m.atom_count}")
    
    # Compact
    compacted = pool.compact(force=True)
    
    # Debug info after compact  
    print(f"After compact: {len(pool.residue_metadata)} metadata entries, compacted={compacted}")
    
    assert compacted == 1
    assert len(pool.residue_metadata) == 2  # Check metadata instead
    assert pool.state.activeResidueCount == 2
    assert pool.stats['compactions'] == 1
    assert pool.get_fragmentation() == 0.0


def test_active_pool_energy_conservation():
    """Test that energy calculations respect active flags."""
    # Create forcefield
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.ljSigma = [0.315, 0.0]
    ff.ljEps = [0.636, 0.0]
    
    pool = ActivePool(box=[5.0, 5.0, 5.0], cutoff=1.2, forcefield=ff)
    
    # Insert interacting molecules
    res0 = pool.insert_molecule(create_water_molecule(2.0, 2.0, 2.0))
    res1 = pool.insert_molecule(create_water_molecule(2.5, 2.0, 2.0))  # Close to first
    
    # Calculate initial energy (simplified for demonstration)
    e_initial = total_energy_active(pool)
    assert e_initial != 0.0  # Should have interaction
    
    # Delete one molecule
    pool.delete_residue(res0)
    
    # Energy should drop to zero (only one active molecule)
    e_after_delete = total_energy_active(pool)
    assert e_after_delete == 0.0  # Single molecule has no self-interaction
    
    # Compact and verify energy remains the same
    pool.compact()
    e_after_compact = total_energy_active(pool)
    assert abs(e_after_compact - e_after_delete) < 1e-10


def test_batched_operations():
    """Test batched insert/delete operations."""
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.ljSigma = [0.315, 0.0]
    ff.ljEps = [0.636, 0.0]
    
    pool = BatchedActivePool(
        box=[10.0, 10.0, 10.0], 
        cutoff=2.5, 
        forcefield=ff,
        batch_size=10
    )
    
    # Queue multiple operations
    for i in range(15):
        x = random.uniform(1.0, 9.0)
        y = random.uniform(1.0, 9.0)
        z = random.uniform(1.0, 9.0)
        pool.queue_insert(create_water_molecule(x, y, z))
    
    # Should not execute yet (batch size is 10, but we have 15 inserts)
    assert pool.state.activeResidueCount == 0
    
    # Check if should flush
    assert pool.should_flush()
    
    # Execute batch
    n_ins, n_del = pool.execute_batch()
    assert n_ins == 15
    assert n_del == 0
    assert pool.state.activeResidueCount == 15
    
    # Queue some deletes
    for i in [0, 2, 4, 6, 8]:
        pool.queue_delete(i)
    
    # Execute deletes
    n_ins, n_del = pool.execute_batch()
    assert n_ins == 0
    assert n_del == 5
    assert pool.state.activeResidueCount == 10


def test_fragmentation_management():
    """Test automatic compaction based on fragmentation threshold."""
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.ljSigma = [0.315, 0.0]
    ff.ljEps = [0.636, 0.0]
    
    pool = BatchedActivePool(
        box=[10.0, 10.0, 10.0],
        cutoff=2.5,
        forcefield=ff,
        batch_size=5
    )
    
    # Insert 10 molecules
    for i in range(10):
        pool.queue_insert(create_water_molecule(i*0.5, 2.5, 2.5))
    pool.execute_batch()
    
    initial_count = pool.state.activeResidueCount
    assert initial_count == 10
    
    # Delete 60% to trigger auto-compact (threshold is 50%)
    for i in range(6):
        pool.queue_delete(i)
    
    # Execute should trigger auto-compact
    n_ins, n_del = pool.execute_batch()
    assert n_del == 6
    
    # After compaction, should have 4 residues
    assert len(pool.residue_metadata) == 4  # Check metadata instead
    assert pool.state.activeResidueCount == 4
    assert pool.stats['compactions'] >= 1


def test_gpu_ready_capacity_limits():
    """Test that capacity limits are properly set for GPU pre-allocation."""
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.ljSigma = [0.315, 0.0]
    ff.ljEps = [0.636, 0.0]
    
    # Create pool with specific capacity
    max_atoms = 100000
    max_residues = 30000
    
    pool = ActivePool(
        box=[20.0, 20.0, 20.0],
        cutoff=2.5,
        forcefield=ff,
        max_atoms=max_atoms,
        max_residues=max_residues
    )
    
    # Verify capacity hints are set
    assert pool.max_atoms == max_atoms
    assert pool.max_residues == max_residues
    
    # These values will be used by GPU kernels for pre-allocation


def test_energy_component_separation():
    """Test separate tracking of vdW and electrostatic energies."""
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.ljSigma = [0.315, 0.0]
    ff.ljEps = [0.636, 0.0]
    
    pool = ActivePool(box=[5.0, 5.0, 5.0], cutoff=1.2, forcefield=ff)
    
    # Insert molecules
    pool.insert_molecule(create_water_molecule(2.0, 2.0, 2.0))
    pool.insert_molecule(create_water_molecule(2.5, 2.0, 2.0))
    
    # Calculate energies (simplified)
    vdw, elec, total = calculate_active_energy_components(pool)
    
    # Verify component sum
    assert abs(total - (vdw + elec)) < 1e-10
    
    # For water with these parameters, should have non-zero vdW
    assert vdw != 0.0


def test_concurrent_insert_delete_pattern():
    """Test realistic GCMC pattern with concurrent inserts and deletes."""
    ff = pygcmc.MCForceField()
    ff.numTotalTypes = 2
    ff.ljSigma = [0.315, 0.0]
    ff.ljEps = [0.636, 0.0]
    
    pool = ActivePool(box=[10.0, 10.0, 10.0], cutoff=2.5, forcefield=ff)
    
    # Simulate GCMC steps
    n_steps = 100
    target_molecules = 20
    
    for step in range(n_steps):
        current_count = pool.state.activeResidueCount
        
        # Decide move type based on current state
        if current_count < target_molecules:
            # Bias toward insertion
            if random.random() < 0.7:
                x = random.uniform(1.0, 9.0)
                y = random.uniform(1.0, 9.0)
                z = random.uniform(1.0, 9.0)
                pool.insert_molecule(create_water_molecule(x, y, z))
        else:
            # Bias toward deletion
            if random.random() < 0.7 and current_count > 0:
                active_indices = pool.get_active_residues()
                if active_indices:
                    idx = random.choice(active_indices)
                    pool.delete_residue(idx)
        
        # Periodic compaction
        if step % 20 == 0 and pool.get_fragmentation() > 0.3:
            pool.compact()
    
    # Verify system reached reasonable state
    assert 10 <= pool.state.activeResidueCount <= 30
    assert pool.stats['total_inserts'] > 0
    assert pool.stats['total_deletes'] > 0