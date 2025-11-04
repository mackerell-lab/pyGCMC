"""
Precise testing of PME Total offset relationship with residues

Hypothesis: offset = f(number of residues)
"""

import sys
import os
import subprocess
import textwrap

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pytest
import pygcmc
pygcmc.set_platform_verbose(True)
pygcmc.set_platform_log_level(pygcmc.PlatformLogLevel.DEBUG)
from pygcmc import MCState, MCAtom, MCResidue, MCForceField
pygcmc.set_platform_verbose(True)
pygcmc.set_platform_log_level(pygcmc.PlatformLogLevel.DEBUG)


def compute_pme_with_state_isolation(n_residues, residue_config):
    """Calculate PME energy in isolated subprocess to avoid global state pollution"""
    
    mesh_size = [32, 32, 32]

    # Build Python script string that directly returns results
    script = textwrap.dedent(f"""
        import pygcmc
        from pygcmc import MCState, MCAtom, MCResidue, MCForceField

        box_size = 5.0
        cutoff = 2.0
        alpha = 2.5
        mesh_size = {mesh_size}
        spline_order = 4

        state = MCState()
        state.info.box = [box_size, box_size, box_size]
        state.info.cutoff = cutoff

        ff = MCForceField()
        ff.numTotalTypes = 1
        ff.numMovementTypes = 1
        ff.ljEps = [0.0]
        ff.ljSigma = [0.3]
        state.forcefield = ff

        positions = [[2.0, 2.0, 2.5], [3.0, 2.0, 2.5], [3.0, 3.0, 2.5], [2.0, 3.0, 2.5]]
        charges = [1.0, -1.0, 1.0, -1.0]

        for pos, charge in zip(positions, charges):
            atom = MCAtom()
            atom.x, atom.y, atom.z = pos
            atom.charge = charge
            atom.type = 0
            state.atoms.append(atom)

        state.activeAtomCount = len(state.atoms)

        residue_config = {residue_config}
        for res_info in residue_config:
            res = MCResidue()
            res.active = True
            res.fixed = False
            res.atomStart = res_info['atomStart']
            res.atomCount = res_info['atomCount']
            res.type = 0
            state.residues.append(res)

        state.activeResidueCount = len(state.residues)

        pygcmc.initializePMEParameters(
            cutoff,
            [box_size, box_size, box_size],
            alpha,
            mesh_size,
            spline_order
        )

        pygcmc.computeSystemEnergyPME(state)

        print({{
            'total': state.ewald_energy.get('total', 0.0),
            'real_space': state.ewald_energy.get('real_space', 0.0),
            'reciprocal': state.ewald_energy.get('reciprocal', 0.0),
            'self': state.ewald_energy.get('self', 0.0)
        }})
    """)

    # Run in subprocess
    try:
        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
        build_path = os.path.join(project_root, 'build')
        existing_path = os.environ.get('PYTHONPATH', '')
        combined_pythonpath = build_path if not existing_path else os.pathsep.join([build_path, existing_path])
        result = subprocess.run(
            [sys.executable, '-c', script],
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, 'PYTHONPATH': combined_pythonpath}
        )

        print("SUBPROCESS STDOUT:", result.stdout)
        print("SUBPROCESS STDERR:", result.stderr)
        result.check_returncode()
        return eval(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running subprocess: {e.stderr}")
        raise


def test_residue_offset_pattern():
    """Test that PME Total should not have residue count-related offset (regression test)
    
    Verify that fixed PME total calculates correctly under different residue configurations
    Create new MCState instances for each test case to avoid state pollution
    """
    
    # Test case 1: 1 residue with 4 atoms
    print("Test case 1: 1 residue with 4 atoms")
    residue_config1 = [{'atomStart': 0, 'atomCount': 4}]
    result1 = compute_pme_with_state_isolation(1, residue_config1)
    
    total1 = result1['total']
    sum1 = result1['real_space'] + result1['reciprocal'] + result1['self']
    offset1 = total1 - sum1
    
    print(f"  Total: {total1:.6f}, Sum: {sum1:.6f}, Offset: {offset1:.6f}")
    assert abs(offset1) < 1e-6, f"Offset for 1 residue configuration should be 0, actual: {offset1}"
    
    # Test case 2: 2 residues with 2 atoms each
    print("\nTest case 2: 2 residues with 2 atoms each")
    residue_config2 = [
        {'atomStart': 0, 'atomCount': 2},
        {'atomStart': 2, 'atomCount': 2}
    ]
    result2 = compute_pme_with_state_isolation(2, residue_config2)
    
    total2 = result2['total']
    sum2 = result2['real_space'] + result2['reciprocal'] + result2['self']
    offset2 = total2 - sum2
    
    print(f"  Total: {total2:.6f}, Sum: {sum2:.6f}, Offset: {offset2:.6f}")
    assert abs(offset2) < 1e-6, f"Offset for 2 residue configuration should be 0, actual: {offset2}"
    
    # Test case 3: 4 residues with 1 atom each
    print("\nTest case 3: 4 residues with 1 atom each")
    residue_config3 = [
        {'atomStart': 0, 'atomCount': 1},
        {'atomStart': 1, 'atomCount': 1},
        {'atomStart': 2, 'atomCount': 1},
        {'atomStart': 3, 'atomCount': 1}
    ]
    result3 = compute_pme_with_state_isolation(4, residue_config3)
    
    total3 = result3['total']
    sum3 = result3['real_space'] + result3['reciprocal'] + result3['self']
    offset3 = total3 - sum3
    
    print(f"  Total: {total3:.6f}, Sum: {sum3:.6f}, Offset: {offset3:.6f}")
    assert abs(offset3) < 1e-6, f"Offset for 4 residue configuration should be 0, actual: {offset3}"
    
    # Verify all offsets are 0
    print(f"\nAll tests passed! Offsets: {offset1:.6f}, {offset2:.6f}, {offset3:.6f}")
    assert abs(offset1) < 1e-6 and abs(offset2) < 1e-6 and abs(offset3) < 1e-6, \
        "Fixed PME Total should not have residue count-related offset"


if __name__ == "__main__":
    test_residue_offset_pattern()
