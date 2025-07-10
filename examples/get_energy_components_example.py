#!/usr/bin/env python3
"""
Example demonstrating how to use getTotalEnergyComponents to get 
separate electrostatic and van der Waals energy components.
"""

import pygcmc

# Create a simple two-particle system
state = pygcmc.MCState()
state.info.cutoff = 12.0
state.info.box = [30.0, 30.0, 30.0]

# Create two atoms with opposite charges
atom1 = pygcmc.MCAtom()
atom1.x, atom1.y, atom1.z = 0.0, 0.0, 0.0
atom1.charge = 1.0  # Positive charge
atom1.type = 0

atom2 = pygcmc.MCAtom()
atom2.x, atom2.y, atom2.z = 2.5, 0.0, 0.0  # 2.5 Angstroms away
atom2.charge = -1.0  # Negative charge
atom2.type = 0

state.atoms = [atom1, atom2]
state.activeAtomCount = 2

# Create residues (one atom per residue)
res1 = pygcmc.MCResidue()
res1.atomStart = 0
res1.atomCount = 1
res1.active = True

res2 = pygcmc.MCResidue()
res2.atomStart = 1
res2.atomCount = 1
res2.active = True

state.residues = [res1, res2]
state.activeResidueCount = 2

# Set up Lennard-Jones parameters
state.forcefield.numTotalTypes = 1
state.forcefield.ljSigma = [2.0]  # sigma = 2.0 Angstroms
state.forcefield.ljEps = [1.0]    # epsilon = 1.0 kJ/mol

# Calculate the system energy
pygcmc.computeSystemEnergy(state)

# Get the energy components separately
elec_energy, vdw_energy = pygcmc.getTotalEnergyComponents(state)

print(f"System Energy Components:")
print(f"  Electrostatic energy: {elec_energy:.4f} kJ/mol")
print(f"  Van der Waals energy: {vdw_energy:.4f} kJ/mol")
print(f"  Total energy: {elec_energy + vdw_energy:.4f} kJ/mol")

# You can also examine individual residue energies
print(f"\nIndividual Residue Energies (before division by 2):")
for i, res in enumerate(state.residues):
    if res.active:
        print(f"  Residue {i}: elec={res.energy_elec:.4f}, vdw={res.energy_vdw:.4f}")

# Note: The residue energies are doubled because each pairwise interaction
# is added to both residues. The getTotalEnergyComponents function 
# automatically divides by 2 to get the correct system total.