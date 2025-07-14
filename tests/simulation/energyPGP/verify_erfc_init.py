# tests/simulation/energyPGP/verify_erfc_init.py
"""
Verify the erfc initialization issue.
"""

import sys
sys.path.insert(0, '/home/zhaomt/gcmc/test107/pygcmc_dev/build')
import pygcmc
from . import pgp_wrapper
from .pgp_wrapper import setPMEParameters, initializePMEParameters, computeSystemEnergyPME
import math
from pygcmc import MCState, MCAtom, MCResidue
from pygcmc import MCForceField

def verify_erfc_init():
    """Verify erfc initialization."""
    print("\n=== Verifying erfc Initialization ===")
    
    # Check if the issue is in the initialization order
    for test_num in range(1, 3):
        print(f"\n--- Test {test_num} ---")
        
        # Create minimal system
        state = MCState()
        state.info.box = [10.0, 10.0, 10.0]
        state.info.cutoff = 1.2
        
        ff = MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff
        
        atoms = []
        atom1 = MCAtom()
        atom1.x = 5.0
        atom1.y = 5.0
        atom1.z = 5.0
        atom1.charge = 1.0
        atom1.type = 0
        atoms.append(atom1)
        
        atom2 = MCAtom()
        atom2.x = 6.0
        atom2.y = 5.0
        atom2.z = 5.0
        atom2.charge = -1.0
        atom2.type = 0
        atoms.append(atom2)
        
        state.atoms = atoms
        state.activeAtomCount = 2
        
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
        
        # Try different initialization orders
        alpha = 2.5
        mesh_size = [32, 32, 32]
        
        if test_num == 1:
            print("Order: setPMEParameters -> initializePMEParameters")
            setPMEParameters(alpha, mesh_size, 4, 1e-6)
            initializePMEParameters(state.info.cutoff, state.info.box, alpha)
        else:
            print("Order: initializePMEParameters only (auto mode)")
            initializePMEParameters(state.info.cutoff)  # alpha=0 for auto
        
        computeSystemEnergyPME(state)
        
        print(f"Total energy: {state.ewald_energy.get('total'):.2f} kJ/mol")
        print(f"Self energy: {state.ewald_energy.get('self', 0.0):.2f} kJ/mol")
        
        # The self energy is correct, so alpha is being set
        # But the residue energies are wrong, so erfcApprox is broken
    
    print("\n\nCONCLUSION:")
    print("The issue is that erfcApprox is returning a constant value.")
    print("This suggests:")
    print("1. The erfc table is not being properly initialized")
    print("2. OR the table lookup is accessing the wrong index")
    print("3. OR there's a default/uninitialized value being returned")
    
    # Let's check what the expected energy should be
    print("\n\nExpected energy calculation:")
    kC = 138.935456
    r = 1.0
    alpha = 2.5
    erfc_val = math.erfc(alpha * r)
    
    print(f"For r=1.0 nm, alpha=2.5:")
    print(f"  erfc(αr) = {erfc_val:.6f}")
    print(f"  Direct Coulomb: -kC/r = {-kC:.2f} kJ/mol")
    print(f"  Real-space only: -kC*erfc(αr)/r = {-kC*erfc_val:.6f} kJ/mol")
    
    print("\nThe fact that we're getting -19303 kJ/mol per residue suggests:")
    print(f"  erfcApprox is returning {-19303 * 2 / (-kC):.6f} instead of {erfc_val:.6f}")
    print(f"  That's a factor of {(-19303 * 2 / (-kC)) / erfc_val:.1f}")

if __name__ == "__main__":
    verify_erfc_init()
