# tests/simulation/energyPGP/test_wrapper_working.py
"""
Test if wrapper is properly loaded and working.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
from pygcmc import setPMEParameters, initializePMEParameters, computeSystemEnergyPGP
from pygcmc import setPGPParameters, precomputeGridPotential


def test_wrapper_loaded():
    """Test if our wrapper functions are loaded."""
    print("\n=== Testing Wrapper Loading ===")
    print(f"computeSystemEnergyPGP type: {type(pygcmc.computeSystemEnergyPGP)}")
    print(f"Is it a function? {hasattr(pygcmc.computeSystemEnergyPGP, '__call__')}")
    print(f"Has __name__? {hasattr(pygcmc.computeSystemEnergyPGP, '__name__')}")
    
    if hasattr(pygcmc.computeSystemEnergyPGP, '__name__'):
        print(f"Function name: {pygcmc.computeSystemEnergyPGP.__name__}")
    
    print(f"\n_orig_computeSystemEnergyPGP exists? {hasattr(pygcmc, '_orig_computeSystemEnergyPGP')}")
    

def test_lj_calculation():
    """Test if LJ energy is calculated without cap."""
    print("\n=== Testing LJ Calculation ===")
    
    # Create state with very close atoms
    state = MCState()
    state.info.box = [10.0, 10.0, 10.0]
    state.info.cutoff = 2.0
    
    ff = MCForceField()
    ff.numTotalTypes = 1
    ff.numMovementTypes = 1
    ff.ljEps = [1.0]
    ff.ljSigma = [0.34]
    state.forcefield = ff
    
    # Two atoms very close
    distance = 0.068  # 0.2 * sigma
    atoms = []
    
    atom1 = MCAtom()
    atom1.x = 5.0
    atom1.y = 5.0
    atom1.z = 5.0
    atom1.charge = 0.0
    atom1.type = 0
    atoms.append(atom1)
    
    atom2 = MCAtom()
    atom2.x = 5.0 + distance
    atom2.y = 5.0
    atom2.z = 5.0
    atom2.charge = 0.0
    atom2.type = 0
    atoms.append(atom2)
    
    state.atoms = atoms
    state.activeAtomCount = 2
    
    # Two residues
    residues = []
    for i in range(2):
        res = MCResidue()
        res.active = True
        res.fixed = False
        res.atomStart = i
        res.atomCount = 1
        res.type = 0
        residues.append(res)
    
    state.residues = residues
    state.activeResidueCount = 2
    
    # Initialize and compute
    alpha = 2.5
    setPMEParameters(alpha, [32, 32, 32], 4, 1e-6)
    initializePMEParameters(state.info.cutoff, state.info.box, alpha)
    setPGPParameters(alpha, [32, 32, 32], state.info.cutoff, [32, 32, 32], 4, 1e-6)
    precomputeGridPotential(state, fixed_only=False)
    
    print(f"Before compute: res[0].energy_vdw = {state.residues[0].energy_vdw}")
    
    computeSystemEnergyPGP(state)
    
    print(f"After compute: res[0].energy_vdw = {state.residues[0].energy_vdw}")
    print(f"After compute: res[1].energy_vdw = {state.residues[1].energy_vdw}")
    
    # Calculate expected
    sigma = 0.34
    epsilon = 1.0
    r_ratio = sigma / distance
    expected = 4 * epsilon * (r_ratio**12 - r_ratio**6)
    
    print(f"\nExpected LJ energy: {expected:.2e}")
    print(f"Total VDW from residues: {state.residues[0].energy_vdw + state.residues[1].energy_vdw}")
    print(f"Per pair (÷2): {(state.residues[0].energy_vdw + state.residues[1].energy_vdw)/2}")
    

if __name__ == "__main__":
    test_wrapper_loaded()
    test_lj_calculation()