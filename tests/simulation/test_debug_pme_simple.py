import pygcmc
from pygcmc import MCState, MCAtom, MCResidue, MCForceField

def test_create_system():
    """Debug test for creating system"""
    state = MCState()
    print(f"Initial state: atoms={state.atoms}, residues={state.residues}")
    print(f"Initial counts: activeAtomCount={state.activeAtomCount}, activeResidueCount={state.activeResidueCount}")
    
    # Try to initialize arrays
    state.atoms = []
    state.residues = []
    print(f"After init: atoms={state.atoms}, residues={state.residues}")
    
    # Create atom
    atom = MCAtom()
    atom.x = 1.5
    atom.y = 1.5
    atom.z = 1.5
    atom.charge = 1.0
    atom.type = 0
    
    # Try different methods to add atom
    print(f"\nAtom created: x={atom.x}, charge={atom.charge}")
    
    # Method 1: append
    try:
        state.atoms.append(atom)
        print(f"Method 1 (append) success: len(atoms)={len(state.atoms)}")
    except Exception as e:
        print(f"Method 1 (append) failed: {e}")
    
    # Method 2: use addAtom
    try:
        idx = state.addAtom(atom)
        print(f"Method 2 (addAtom) success: returned index={idx}, activeAtomCount={state.activeAtomCount}")
    except Exception as e:
        print(f"Method 2 (addAtom) failed: {e}")
    
    print(f"\nFinal state: len(atoms)={len(state.atoms)}, activeAtomCount={state.activeAtomCount}")
    print(f"Atoms array: {state.atoms}")

if __name__ == "__main__":
    test_create_system()