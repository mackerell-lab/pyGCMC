# PyGCMC Energy Functions Analysis for GCMC

## Available PyGCMC Energy Functions

### System Energy Functions
1. `computeSystemEnergy()` - Full system energy, no cutoff
2. `computeSystemEnergyCutoff()` - System energy with distance cutoff
3. `computeSystemEnergyPBC()` - System energy with periodic boundaries
4. `computeSystemEnergyPBCCutoff()` - PBC with cutoff
5. `computeSystemVdwEnergyCutoff()` - Only VDW component
6. `computeSystemEnergyEwald()` - Ewald summation for long-range
7. `computeSystemEnergyPME()` - PME for long-range electrostatics
8. `computeSystemEnergyPGP()` - Planar Gaussian Potential

### Movement Energy Functions
1. `computeMovementEnergy()` - Energy for movement residues only
2. `computeMovementEnergyCutoff()` - Movement energy with cutoff
3. `computeMovementEnergyEwald()` - Movement with Ewald
4. `computeMovementEnergyPME()` - Movement with PME
5. `computeMovementEnergyPGP()` - Movement with PGP

## GCMC Requirements

### For Insertion/Deletion:
1. **Calculate energy of single molecule** - Need to compute interaction energy between:
   - New molecule ↔ Protein atoms
   - New molecule ↔ Existing guest molecules
   - New molecule ↔ Water/ions

2. **Energy difference calculation**:
   - ΔE = E_new - E_old
   - For insertion: E_old is system without molecule
   - For deletion: E_new is system without molecule

### For Translation/Rotation:
1. **Single molecule energy before/after move**
2. **Only interactions with other molecules change**

## Analysis: Are Current Functions Sufficient?

### ✅ What we have:
1. **Movement energy functions** - Perfect for GCMC!
   - Can calculate energy of just the guest molecules
   - Avoids recalculating protein-protein interactions
   - `computeMovementEnergyCutoff()` is ideal

2. **System energy** - For validation and total energy

3. **Energy storage** - Each residue stores:
   - `residue.energy_vdw`
   - `residue.energy_elec`

### ❓ What might be missing:

1. **Single residue energy calculation**:
   - Current functions calculate all movement residues
   - GCMC needs energy of ONE specific molecule
   - Workaround: Set only target residue as movement type temporarily

2. **Partial energy updates**:
   - After move, only need to update moved molecule's energy
   - Current: Must recalculate all movement residues

3. **Cavity detection**:
   - Need to check if insertion point is valid
   - Can use distance calculations

## Recommended Approach

### For GCMC Insert Test:
```python
def calculate_insertion_energy(system, new_molecule_atoms, position):
    # Create trial system with new molecule
    trial_system = system.copy()
    
    # Add molecule at position
    residue = add_molecule_to_system(trial_system, new_molecule_atoms, position)
    
    # Option 1: Use movement energy
    # Set new residue as only movement residue
    trial_system.residue.type = 0  # movement
    pygcmc.computeMovementEnergyCutoff(trial_system)
    
    # Get energy from residue
    insertion_energy = residue.energy_vdw + residue.energy_elec
    
    return insertion_energy
```

### For GCMC with existing molecules:
```python
def calculate_delta_energy(system_before, system_after):
    # Calculate total energies
    pygcmc.computeSystemEnergyCutoff(system_before)
    pygcmc.computeSystemEnergyCutoff(system_after)
    
    # Sum energies (divide by 2 for double counting)
    E_before = sum(res.energy_vdw + res.energy_elec for res in system_before.residues) / 2
    E_after = sum(res.energy_vdw + res.energy_elec for res in system_after.residues) / 2
    
    return E_after - E_before
```

## Conclusion

✅ **PyGCMC has sufficient energy functions for GCMC implementation!**

The movement energy functions are particularly well-suited for GCMC because they:
1. Calculate only guest-host and guest-guest interactions
2. Store results per residue
3. Support various methods (cutoff, PME, Ewald)

We can proceed with the Python prototype using these functions!