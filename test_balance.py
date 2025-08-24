import pygcmc
import numpy as np

# Setup
state = pygcmc.MCState()
state.info.box = np.array([3.0, 3.0, 3.0])

ff = pygcmc.MCForceField()
ff.numTotalTypes = 1
ff.numMovementTypes = 1
ff.ljEps = [0.5]
ff.ljSigma = [0.3]
state.forcefield = ff

params = pygcmc.movement.MovementParams()
params.temperature = 300.0
params.chemicalPotential = -15.0
params.seed = 42

mover = pygcmc.movement.MovementModule()
mover.setParams(params)

# Track paired operations
paired_insertions = 0
paired_deletions = 0

for i in range(100):
    result_ins = mover.attemptInsertion(state)
    if result_ins.accepted:
        paired_insertions += 1
        # Immediately try deletion - should prefer the just-inserted residue
        result_del = mover.attemptDeletion(state)
        if result_del.accepted:
            paired_deletions += 1
            print(f"Pair {i}: Insert at idx {result_ins.residueIndex}, Delete at idx {result_del.residueIndex}")
        else:
            print(f"Pair {i}: Insert at idx {result_ins.residueIndex}, Delete REJECTED")

print(f"\nPaired insertions: {paired_insertions}")
print(f"Paired deletions: {paired_deletions}")
if paired_insertions > 0:
    print(f"Ratio: {paired_insertions/max(1, paired_deletions):.2f}")
