"""Debug deletion probability issue"""

import os
import pygcmc

# Create system exactly as in test
state = pygcmc.MCState()
state.info.box = (30.0, 30.0, 30.0)
state.info.setTemperature(300.0)

# Zero interactions for ideal gas
ff = pygcmc.MCForceField()
ff.numTotalTypes = 1
ff.numMovementTypes = 1
ff.ljSigma = [0.315]  # nm
ff.ljEps = [0.0]      # Zero interaction - ideal gas
state.forcefield = ff

# Single atom template
template = pygcmc.movement.FragmentTemplate()
template.typeId = 0
template.atoms = [pygcmc.MCAtom()]
template.atoms[0].type = 0
template.atoms[0].charge = 0.0

reservoir = pygcmc.movement.FragmentReservoir()
reservoir.addTemplate(template)

engine = pygcmc.GCMCEngine()
engine.initialize(state, reservoir)
engine.setTemperature(300.0)
engine.setSeed(12345)

box_volume = 30.0 * 30.0 * 30.0  # = 27000 nm^3

acceptance = pygcmc.GCMCAcceptance()
acceptance.setTemperature(300.0)
acceptance.setVolume(box_volume)
acceptance.setActivity(0, 0.001)
acceptance.setSeed(12346)
engine.setAcceptanceCalculator(acceptance)

# Enable probability storage
os.environ['GCMC_STORE_PROB'] = '1'

# Insert exactly 10 molecules
count = 0
attempts = 0
while count < 10 and attempts < 1000:
    result = engine.attemptInsertion(0)
    if result.accepted:
        count += 1
    attempts += 1

print(f"Inserted {count} molecules after {attempts} attempts")
print(f"Current N = {reservoir.getActiveCount(0)}")

# Now try deletion at N=10
del_result = engine.attemptDeletion(0)
print(f"\nDeletion at N=10:")
print(f"  acceptanceProbability = {del_result.acceptanceProbability}")
print(f"  accepted = {del_result.accepted}")

# Expected value
zV = 0.001 * 27000.0  # = 27.0
expected = min(1.0, 10 / zV)  # = 0.370370
print(f"  expected = {expected:.10f}")
print(f"  ratio = {del_result.acceptanceProbability / expected:.10f}")

# Check what acceptance calculator says directly
direct_prob = acceptance.calculateDeletionProbability(0, 10, 0.0, 1.0)
print(f"\nDirect from acceptance calculator:")
print(f"  prob = {direct_prob:.10f}")