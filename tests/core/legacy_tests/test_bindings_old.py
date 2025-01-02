# tests/core/test_bindings.py

import pyGCMC_bindings

# Create a System instance
sys = pyGCMC_bindings.System(epsilon=1.0, sigma=3.355)

# Create Particles
h1 = pyGCMC_bindings.Particle(serial=1, name="H1", residue="HOH", sequence=1, x=0.0, y=0.0, z=0.0, charge=0.5, type="H", nameTop="H")
o1 = pyGCMC_bindings.Particle(serial=2, name="O1", residue="HOH", sequence=1, x=1.0, y=0.0, z=0.0, charge=-1.0, type="O", nameTop="O")
h2 = pyGCMC_bindings.Particle(serial=3, name="H2", residue="HOH", sequence=1, x=1.0, y=1.0, z=0.0, charge=0.5, type="H", nameTop="H")

# Create a Residue
water = pyGCMC_bindings.Residue(name="HOH", sequence_number=1, chain_id='A')
water.atoms = [h1, o1, h2]

# Add Residue to System
sys.add_residue(water)

# Print system state
kinetic, potential = sys.get_system_state()
print(f"Kinetic Energy: {kinetic}")
print(f"Potential Energy: {potential}")
print(f"Total Energy: {kinetic + potential}")

# Access residues and their particles
for i in range(sys.get_residue_count()):
    residue = sys.get_residue(i)
    print(f"Residue {i}: {residue.name}, Sequence: {residue.sequence_number}, Chain: {residue.chain_id}")
    for atom in residue.atoms:
        print(f"  Particle {atom.serial}: {atom.name}, Position: ({atom.x}, {atom.y}, {atom.z})")

# Update positions and velocities
sys.update_positions(0.001)
sys.update_velocities(0.001)

# Print updated system state
kinetic, potential = sys.get_system_state()
print(f"After update - Kinetic Energy: {kinetic}")
print(f"After update - Potential Energy: {potential}")
print(f"After update - Total Energy: {kinetic + potential}")
