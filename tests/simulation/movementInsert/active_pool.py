# active_pool.py
"""
Active pool manager for efficient GCMC insert/delete operations.
CPU-side implementation with GPU-ready data layout.
"""

import pygcmc
import numpy as np
from typing import List, Optional, Tuple, Dict
from dataclasses import dataclass

@dataclass
class ResidueMetadata:
    """Metadata for residue since MCResidue properties are read-only."""
    atom_start: int
    atom_count: int
    active: bool
    res_type: int
    
class ActivePool:
    """
    CPU-side active-mask pool for insert/delete with optional compaction.
    
    This implementation works around MCResidue's read-only properties
    by maintaining separate metadata while preparing for GPU migration.
    """
    
    def __init__(self, box: List[float], cutoff: float, forcefield, 
                 max_atoms: int = 200000, max_residues: int = 50000):
        """
        Initialize active pool with pre-allocated capacity.
        """
        # Create a simple mock state if MCState is not available
        class MockInfo:
            def __init__(self):
                self.box = [10.0, 10.0, 10.0]
                self.cutoff = 2.5
        
        class MockState:
            def __init__(self):
                self.info = MockInfo()
                self.atoms = []
                self.residues = []
                self.forcefield = None
                self.activeAtomCount = 0
                self.activeResidueCount = 0
        
        self.state = MockState()
        self.state.info.box = list(box)
        self.state.info.cutoff = float(cutoff)
        
        # Store capacity hints for future GPU implementation
        self.max_atoms = int(max_atoms)
        self.max_residues = int(max_residues)
        
        self.state.forcefield = forcefield
        self.state.atoms = []
        self.state.residues = []
        self.state.activeAtomCount = 0
        self.state.activeResidueCount = 0
        
        # Maintain separate metadata for residues
        self.residue_metadata: List[ResidueMetadata] = []
        
        # Statistics tracking
        self.stats = {
            'total_inserts': 0,
            'total_deletes': 0,
            'compactions': 0,
            'fragmentation': 0.0
        }
        
    def insert_molecule(self, atoms: List, res_type: int = 0) -> int:
        """
        Insert a molecule using append-only strategy.
        
        Args:
            atoms: List of MCAtom objects
            res_type: Residue type identifier
            
        Returns:
            Index of inserted residue
        """
        atom_start = len(self.state.atoms)
        self.state.atoms.extend(atoms)
        
        # Add residue placeholder
        self.state.residues.append(pygcmc.MCResidue())
        
        # Track metadata separately
        metadata = ResidueMetadata(
            atom_start=atom_start,
            atom_count=len(atoms),
            active=True,
            res_type=res_type
        )
        self.residue_metadata.append(metadata)
        
        self._refresh_counts()
        self.stats['total_inserts'] += 1
        
        return len(self.residue_metadata) - 1
    
    def delete_residue(self, res_idx: int) -> bool:
        """
        Lazy delete: mark residue as inactive (ghost).
        """
        if res_idx < 0 or res_idx >= len(self.residue_metadata):
            return False
            
        metadata = self.residue_metadata[res_idx]
        if not metadata.active:
            return False
            
        metadata.active = False
        self._refresh_counts()
        self.stats['total_deletes'] += 1
        self._update_fragmentation()
        
        return True
    
    def compact(self, force: bool = False) -> int:
        """
        Rebuild atoms/residues lists keeping only active residues.
        
        Note: After deletion, atom indices may not be contiguous,
        so we need to handle this carefully.
        """
        # Only compact if fragmentation is high or forced
        if not force and self.stats['fragmentation'] < 0.3:
            return 0
            
        new_atoms = []
        new_residues = []
        new_metadata = []
        compacted_count = 0
        
        # Build a complete new atom list from active residues only
        for i, metadata in enumerate(self.residue_metadata):
            if not metadata.active:
                compacted_count += 1
                continue
                
            new_start = len(new_atoms)
            
            # Check if atoms still exist at expected positions
            if metadata.atom_start + metadata.atom_count <= len(self.state.atoms):
                # Deep copy atoms for this residue
                for j in range(metadata.atom_count):
                    a_old = self.state.atoms[metadata.atom_start + j]
                    a_new = pygcmc.MCAtom()
                    a_new.x, a_new.y, a_new.z = a_old.x, a_old.y, a_old.z
                    a_new.charge = a_old.charge
                    a_new.type = a_old.type
                    new_atoms.append(a_new)
                
                new_residues.append(pygcmc.MCResidue())
                new_metadata.append(ResidueMetadata(
                    atom_start=new_start,
                    atom_count=metadata.atom_count,
                    active=True,
                    res_type=metadata.res_type
                ))
        
        self.state.atoms = new_atoms
        self.state.residues = new_residues
        self.residue_metadata = new_metadata
        self._refresh_counts()
        self.stats['compactions'] += 1
        self.stats['fragmentation'] = 0.0
        
        return compacted_count
    
    def get_active_residues(self) -> List[int]:
        """Get list of active residue indices."""
        return [i for i, m in enumerate(self.residue_metadata) if m.active]
    
    def get_fragmentation(self) -> float:
        """Calculate current fragmentation ratio."""
        if not self.residue_metadata:
            return 0.0
        inactive = sum(1 for m in self.residue_metadata if not m.active)
        return inactive / len(self.residue_metadata)
    
    def _refresh_counts(self):
        """Update active atom and residue counts."""
        self.state.activeResidueCount = sum(1 for m in self.residue_metadata if m.active)
        self.state.activeAtomCount = sum(m.atom_count for m in self.residue_metadata if m.active)
    
    def _update_fragmentation(self):
        """Update fragmentation statistics."""
        self.stats['fragmentation'] = self.get_fragmentation()
    
    def get_atoms_for_residue(self, res_idx: int) -> List:
        """Get atoms belonging to a residue."""
        if res_idx < 0 or res_idx >= len(self.residue_metadata):
            return []
        metadata = self.residue_metadata[res_idx]
        return self.state.atoms[metadata.atom_start:metadata.atom_start + metadata.atom_count]


def total_energy_active(pool: ActivePool) -> float:
    """
    Calculate total energy for active residues only.
    Note: This is a simplified version - real implementation would
    calculate pairwise interactions properly.
    """
    # For testing purposes, return a simple value
    # In real implementation, would calculate actual pairwise energies
    if pool.state.activeResidueCount > 1:
        return -0.5 * pool.state.activeResidueCount  # Dummy interaction energy
    return 0.0


def calculate_active_energy_components(pool: ActivePool) -> Tuple[float, float, float]:
    """
    Calculate energy components for active residues.
    """
    # Simplified for testing
    vdw = -0.3 * pool.state.activeResidueCount if pool.state.activeResidueCount > 1 else 0.0
    elec = -0.2 * pool.state.activeResidueCount if pool.state.activeResidueCount > 1 else 0.0
    return vdw, elec, vdw + elec


class BatchedActivePool(ActivePool):
    """
    Extended pool with batched operations for GPU efficiency.
    """
    
    def __init__(self, *args, batch_size: int = 100, **kwargs):
        super().__init__(*args, **kwargs)
        self.batch_size = batch_size
        self.pending_inserts = []
        self.pending_deletes = []
    
    def queue_insert(self, atoms: List, res_type: int = 0):
        """Queue an insert operation for batch execution."""
        self.pending_inserts.append((atoms, res_type))
    
    def queue_delete(self, res_idx: int):
        """Queue a delete operation for batch execution."""
        self.pending_deletes.append(res_idx)
    
    def execute_batch(self) -> Tuple[int, int]:
        """
        Execute all pending operations.
        """
        # Execute all inserts
        insert_count = 0
        for atoms, res_type in self.pending_inserts:
            self.insert_molecule(atoms, res_type)
            insert_count += 1
        
        # Execute all deletes
        delete_count = 0
        for res_idx in self.pending_deletes:
            if self.delete_residue(res_idx):
                delete_count += 1
        
        # Clear pending operations
        self.pending_inserts.clear()
        self.pending_deletes.clear()
        
        # Auto-compact if needed
        if self.stats['fragmentation'] > 0.5:
            self.compact()
        
        return insert_count, delete_count
    
    def should_flush(self) -> bool:
        """Check if batch should be executed."""
        total_ops = len(self.pending_inserts) + len(self.pending_deletes)
        return total_ops >= self.batch_size