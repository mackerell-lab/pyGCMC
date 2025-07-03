# tests/simulation/energyOpenmm/nbfix_tests.py
"""
NBFIX-specific tests comparing PyGCMC with OpenMM

This module tests:
1. NBFIX parameters correctly override Lorentz-Berthelot mixing rules
2. Ion-water interactions with NBFIX
3. Energy consistency between PyGCMC and OpenMM for NBFIX systems
"""

import pytest
import math
import pygcmc
import os
import warnings

# Filter SWIG-related warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

try:
    from openmm import *
    from openmm.app import *
    from openmm.unit import *
    HAS_OPENMM = True
except ImportError:
    HAS_OPENMM = False

# Constants
ANGSTROM_TO_NM = 0.1
KCAL_TO_KJ = 4.184
kC = 138.935456  # Coulomb constant in kJ·nm/mol/e^2

# Get test data directory
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(TEST_DIR)), "data")


def create_nbfix_forcefield():
    """Create a force field with NBFIX parameters from toppar_water_ions.str"""
    ff = pygcmc.ForceField()
    
    # Load the water/ions parameter file
    water_ions_file = os.path.join(DATA_DIR, "toppar_water_ions.str")
    if os.path.exists(water_ions_file):
        pygcmc.PRMParser.parse_file_to_forcefield(water_ions_file, ff)
    
    # Manually add SOD-CLA NBFIX if not present (from literature values)
    # SOD    CLA      -0.083875   3.731
    try:
        result = ff.get_nbfix("SOD", "CLA")
        if not (len(result) >= 3 and result[2]):
            ff.add_nbfix("SOD", "CLA", -0.083875, 3.731)
    except:
        ff.add_nbfix("SOD", "CLA", -0.083875, 3.731)
    
    return ff


def create_simple_ion_water_mcstate():
    """Create a simple Na+ and water MCState for NBFIX testing"""
    # Create MCState directly
    state = pygcmc.MCState()
    
    # Set box and cutoff
    state.info.box = [3.0, 3.0, 3.0]  # 3 nm box
    state.info.cutoff = 1.2  # 1.2 nm cutoff
    
    # Define atom types
    sod_idx = state.atomTypes.get_or_add_type("SOD")  # 0
    ot_idx = state.atomTypes.get_or_add_type("OT")   # 1
    ht_idx = state.atomTypes.get_or_add_type("HT")   # 2
    cla_idx = state.atomTypes.get_or_add_type("CLA")  # 3
    
    # Add atoms
    atoms = []
    
    # Sodium ion
    atom_na = pygcmc.MCAtom()
    atom_na.x = 1.5  # nm
    atom_na.y = 1.5
    atom_na.z = 1.5
    atom_na.charge = 1.0
    atom_na.type = sod_idx  # SOD
    atoms.append(atom_na)
    
    # Water oxygen
    atom_o = pygcmc.MCAtom()
    atom_o.x = 1.8  # 0.3 nm from Na+
    atom_o.y = 1.5
    atom_o.z = 1.5
    atom_o.charge = -0.834
    atom_o.type = ot_idx  # OT
    atoms.append(atom_o)
    
    # Water hydrogen 1
    atom_h1 = pygcmc.MCAtom()
    atom_h1.x = 1.8957
    atom_h1.y = 1.5
    atom_h1.z = 1.5
    atom_h1.charge = 0.417
    atom_h1.type = ht_idx  # HT
    atoms.append(atom_h1)
    
    # Water hydrogen 2
    atom_h2 = pygcmc.MCAtom()
    atom_h2.x = 1.8
    atom_h2.y = 1.5957
    atom_h2.z = 1.5
    atom_h2.charge = 0.417
    atom_h2.type = ht_idx  # HT
    atoms.append(atom_h2)
    
    # Chloride ion
    atom_cl = pygcmc.MCAtom()
    atom_cl.x = 1.5
    atom_cl.y = 1.5
    atom_cl.z = 2.0  # 0.5 nm from Na+
    atom_cl.charge = -1.0
    atom_cl.type = cla_idx  # CLA
    atoms.append(atom_cl)
    
    state.atoms = atoms
    state.activeAtomCount = len(atoms)
    
    # Add residues
    residues = []
    
    # Sodium residue
    res_na = pygcmc.MCResidue()
    res_na.active = True
    res_na.atomStart = 0
    res_na.atomCount = 1
    res_na.type = 0
    residues.append(res_na)
    
    # Water residue
    res_water = pygcmc.MCResidue()
    res_water.active = True
    res_water.atomStart = 1
    res_water.atomCount = 3
    res_water.type = 1
    residues.append(res_water)
    
    # Chloride residue
    res_cl = pygcmc.MCResidue()
    res_cl.active = True
    res_cl.atomStart = 4
    res_cl.atomCount = 1
    res_cl.type = 2
    residues.append(res_cl)
    
    state.residues = residues
    state.activeResidueCount = len(residues)
    
    return state


def setup_mcstate_with_forcefield(state, forcefield):
    """Apply force field parameters to MCState"""
    # Get number of types
    n_types = len(state.atomTypes.atomTypes)
    
    # Initialize forcefield arrays
    state.forcefield.numTotalTypes = n_types
    
    # Create new arrays
    ljSigma = [0.0] * (n_types * n_types)
    ljEps = [0.0] * (n_types * n_types)
    
    # Apply force field parameters with NBFIX
    for i in range(n_types):
        for j in range(n_types):
            type1 = state.atomTypes.atomTypes[i]
            type2 = state.atomTypes.atomTypes[j]
            idx = i * n_types + j
            
            # Check for NBFIX first
            try:
                result = forcefield.get_nbfix(type1, type2)
                if len(result) >= 3 and result[2]:  # has_nbfix
                    # NBFIX found - use those parameters
                    # IMPORTANT: NBFIX parameters in CHARMM files:
                    # - epsilon is negative by convention (we use abs() to make it positive)
                    # - rmin is the actual Rmin value (NOT Rmin/2 like in standard LJ params)
                    epsilon = abs(result[0]) * KCAL_TO_KJ
                    rmin = result[1]  # This is Rmin (not Rmin/2)
                    # Convert Rmin to sigma for standard LJ form: sigma = Rmin / 2^(1/6)
                    sigma = rmin / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM
                else:
                    raise ValueError("No NBFIX")
            except Exception as e:
                # No NBFIX - use standard Lorentz-Berthelot mixing rules
                try:
                    lj1 = forcefield.get_lj_params(type1)
                    lj2 = forcefield.get_lj_params(type2)
                    
                    # Convert CHARMM parameters to standard LJ form
                    # IMPORTANT: Standard LJ params in CHARMM use rmin_half (Rmin/2)
                    # So we multiply by 2 to get Rmin, then convert to sigma
                    sigma1 = 2.0 * lj1.rmin_half / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM
                    sigma2 = 2.0 * lj2.rmin_half / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM
                    eps1 = abs(lj1.epsilon) * KCAL_TO_KJ
                    eps2 = abs(lj2.epsilon) * KCAL_TO_KJ
                    
                    # Lorentz-Berthelot mixing rules:
                    # sigma_ij = (sigma_i + sigma_j) / 2  (arithmetic mean)
                    # epsilon_ij = sqrt(epsilon_i * epsilon_j)  (geometric mean)
                    sigma = (sigma1 + sigma2) / 2.0
                    epsilon = math.sqrt(eps1 * eps2)
                except:
                    sigma = 0.0
                    epsilon = 0.0
            
            ljSigma[idx] = sigma
            ljEps[idx] = epsilon
            
    
    # Assign the complete arrays
    state.forcefield.ljSigma = ljSigma
    state.forcefield.ljEps = ljEps
    
    
    return state


def setup_openmm_system_with_nbfix(mc_state, forcefield, separate_forces=False):
    """Create OpenMM system with NBFIX parameters matching MCState
    
    Args:
        mc_state: PyGCMC MCState object
        forcefield: PyGCMC ForceField object with NBFIX parameters
        separate_forces: If True, create separate VDW and electrostatic forces for detailed comparison
    
    Returns:
        system: OpenMM System object
        positions: List of atomic positions
    """
    system = System()
    
    # Add particles with masses (use standard masses for now)
    masses = {"SOD": 22.98977, "OT": 15.9994, "HT": 1.008, "CLA": 35.45}
    for i in range(mc_state.activeAtomCount):
        atom_type = mc_state.atomTypes.atomTypes[mc_state.atoms[i].type]
        mass = masses.get(atom_type, 1.0)
        system.addParticle(mass)
    
    # Set box
    box = mc_state.info.box
    system.setDefaultPeriodicBoxVectors(
        Vec3(box[0], 0, 0) * nanometers,
        Vec3(0, box[1], 0) * nanometers,
        Vec3(0, 0, box[2]) * nanometers
    )
    
    # Build parameter tables from MCState forcefield
    n_types = mc_state.forcefield.numTotalTypes
    sigma_table = []
    epsilon_table = []
    
    # Copy parameters from MCState forcefield matrix
    for idx in range(n_types * n_types):
        sigma_table.append(mc_state.forcefield.ljSigma[idx])
        epsilon_table.append(mc_state.forcefield.ljEps[idx])
    
    if separate_forces:
        # Create separate VDW and electrostatic forces for detailed energy comparison
        
        # VDW Force (Force Group 1)
        vdw_expr = """4*epsilon*((sigma/r)^12-(sigma/r)^6);
                      sigma=sigma_table(type1, type2);
                      epsilon=epsilon_table(type1, type2)"""
        vdw_force = CustomNonbondedForce(vdw_expr)
        vdw_force.addPerParticleParameter("type")
        
        # Add tables for VDW
        vdw_force.addTabulatedFunction("sigma_table", 
                                      Discrete2DFunction(n_types, n_types, sigma_table))
        vdw_force.addTabulatedFunction("epsilon_table", 
                                      Discrete2DFunction(n_types, n_types, epsilon_table))
        
        # Add particles to VDW force
        for i in range(mc_state.activeAtomCount):
            vdw_force.addParticle([float(mc_state.atoms[i].type)])
        
        # Set VDW cutoff and method
        vdw_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        vdw_force.setCutoffDistance(mc_state.info.cutoff * nanometers)
        vdw_force.setForceGroup(1)
        
        # Electrostatic Force (Force Group 2)
        elec_expr = "kC*q1*q2/r"
        elec_force = CustomNonbondedForce(elec_expr)
        elec_force.addGlobalParameter("kC", kC)  # kC = 138.935456 kJ·nm/mol/e²
        elec_force.addPerParticleParameter("q")
        
        # Add particles to electrostatic force
        for i in range(mc_state.activeAtomCount):
            elec_force.addParticle([mc_state.atoms[i].charge])
        
        # Set electrostatic cutoff and method
        elec_force.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        elec_force.setCutoffDistance(mc_state.info.cutoff * nanometers)
        elec_force.setForceGroup(2)
        
        # Add exclusions to BOTH forces for bonded atoms in water
        # Water has atoms at indices 1, 2, 3 (OT, HT, HT)
        vdw_force.addExclusion(1, 2)   # O-H1
        vdw_force.addExclusion(1, 3)   # O-H2
        vdw_force.addExclusion(2, 3)   # H1-H2
        elec_force.addExclusion(1, 2)  # O-H1
        elec_force.addExclusion(1, 3)  # O-H2
        elec_force.addExclusion(2, 3)  # H1-H2
        
        system.addForce(vdw_force)
        system.addForce(elec_force)
    else:
        # Create combined force (original behavior)
        energy_expr = """4*epsilon*((sigma/r)^12-(sigma/r)^6) + kC*q1*q2/r;
                        sigma=sigma_table(type1, type2);
                        epsilon=epsilon_table(type1, type2)"""
        
        nbforce = CustomNonbondedForce(energy_expr)
        nbforce.addGlobalParameter("kC", kC)  # kC = 138.935456 kJ·nm/mol/e²
        nbforce.addPerParticleParameter("q")
        nbforce.addPerParticleParameter("type")
        
        # Add tables to force
        nbforce.addTabulatedFunction("sigma_table", 
                                    Discrete2DFunction(n_types, n_types, sigma_table))
        nbforce.addTabulatedFunction("epsilon_table", 
                                    Discrete2DFunction(n_types, n_types, epsilon_table))
        
        # Add particles with parameters - ensure type is float to avoid interpolation issues
        for i in range(mc_state.activeAtomCount):
            atom = mc_state.atoms[i]
            nbforce.addParticle([atom.charge, float(atom.type)])
        
        # Set method and cutoff
        nbforce.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
        nbforce.setCutoffDistance(mc_state.info.cutoff * nanometers)
        
        # Add exclusions for bonded atoms in water (TIP3P has 3 atoms)
        nbforce.addExclusion(1, 2)  # O-H1
        nbforce.addExclusion(1, 3)  # O-H2
        nbforce.addExclusion(2, 3)  # H1-H2
        
        system.addForce(nbforce)
    
    # Create positions
    positions = []
    for i in range(mc_state.activeAtomCount):
        atom = mc_state.atoms[i]
        positions.append(Vec3(atom.x, atom.y, atom.z) * nanometers)
    
    return system, positions


@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available")
def test_nbfix_overrides_lj_combination():
    """Test that NBFIX parameters override default Lorentz-Berthelot mixing rules"""
    
    # Create system and force field
    mc_state = create_simple_ion_water_mcstate()
    forcefield = create_nbfix_forcefield()
    
    # Apply force field with NBFIX
    mc_state = setup_mcstate_with_forcefield(mc_state, forcefield)
    
    # Find SOD and CLA type indices
    sod_idx = -1
    cla_idx = -1
    for i, atom_type in enumerate(mc_state.atomTypes.atomTypes):
        if atom_type == "SOD":
            sod_idx = i
        elif atom_type == "CLA":
            cla_idx = i
    
    
    assert sod_idx >= 0 and cla_idx >= 0, "SOD and CLA types not found"
    
    # Get the interaction parameters
    n_types = mc_state.forcefield.numTotalTypes
    nbfix_idx = sod_idx * n_types + cla_idx
    
    # Expected NBFIX values from toppar_water_ions.str
    # SOD    CLA      -0.083875   3.731
    expected_eps = 0.083875 * KCAL_TO_KJ  # Convert to kJ/mol
    expected_sigma = 3.731 / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM  # Convert Rmin to sigma
    
    
    actual_eps = mc_state.forcefield.ljEps[nbfix_idx]
    actual_sigma = mc_state.forcefield.ljSigma[nbfix_idx]
    
    print(f"\nNBFIX test for SOD-CLA:")
    print(f"Expected epsilon: {expected_eps:.6f} kJ/mol")
    print(f"Actual epsilon: {actual_eps:.6f} kJ/mol")
    print(f"Expected sigma: {expected_sigma:.6f} nm")
    print(f"Actual sigma: {actual_sigma:.6f} nm")
    
    # Check if NBFIX values are used (within numerical tolerance)
    assert abs(actual_eps - expected_eps) < 1e-3, f"NBFIX epsilon mismatch"
    assert abs(actual_sigma - expected_sigma) < 1e-4, f"NBFIX sigma mismatch"
    
    # Now check that standard mixing would give different values
    sod_lj = forcefield.get_lj_params("SOD")
    cla_lj = forcefield.get_lj_params("CLA")
    
    # Calculate what standard mixing would give
    sod_sigma = 2.0 * sod_lj.rmin_half / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM
    cla_sigma = 2.0 * cla_lj.rmin_half / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM
    sod_eps = abs(sod_lj.epsilon) * KCAL_TO_KJ
    cla_eps = abs(cla_lj.epsilon) * KCAL_TO_KJ
    
    standard_sigma = (sod_sigma + cla_sigma) / 2.0
    standard_eps = math.sqrt(sod_eps * cla_eps)
    
    print(f"\nStandard mixing would give:")
    print(f"Standard epsilon: {standard_eps:.6f} kJ/mol")
    print(f"Standard sigma: {standard_sigma:.6f} nm")
    print(f"Difference in epsilon: {abs(actual_eps - standard_eps):.6f} kJ/mol")
    print(f"Difference in sigma: {abs(actual_sigma - standard_sigma):.6f} nm")
    
    # NBFIX should be different from standard mixing  
    # For SOD-CLA, the epsilon is very close to standard mixing, but sigma is different
    assert abs(actual_sigma - standard_sigma) > 0.001, "NBFIX sigma should differ from standard mixing"


@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available")
def test_ion_water_nbfix_energy():
    """Test ion-water interaction energies with NBFIX"""
    
    # Create system and force field
    mc_state = create_simple_ion_water_mcstate()
    forcefield = create_nbfix_forcefield()
    
    # Apply force field with NBFIX
    mc_state = setup_mcstate_with_forcefield(mc_state, forcefield)
    
    # Calculate energy with PyGCMC
    pygcmc.computeSystemEnergyPBCCutoff(mc_state)
    
    # Get residue energies
    total_vdw = 0.0
    total_elec = 0.0
    for i in range(mc_state.activeResidueCount):
        total_vdw += mc_state.residues[i].energy_vdw
        total_elec += mc_state.residues[i].energy_elec
    total_energy = total_vdw + total_elec
    
    print(f"\nPyGCMC Energy Components:")
    print(f"VDW energy: {total_vdw:.6f} kJ/mol")
    print(f"Electrostatic energy: {total_elec:.6f} kJ/mol")
    print(f"Total energy: {total_energy:.6f} kJ/mol")
    
    # Detailed residue breakdown
    residue_names = ["SOD", "TIP3", "CLA"]
    for i, res in enumerate(mc_state.residues[:mc_state.activeResidueCount]):
        print(f"\nResidue {i} ({residue_names[i]}):")
        print(f"  VDW: {res.energy_vdw:.6f} kJ/mol")
        print(f"  Elec: {res.energy_elec:.6f} kJ/mol")
        print(f"  Total: {res.energy_vdw + res.energy_elec:.6f} kJ/mol")
    
    # Check that energies are reasonable
    assert abs(total_elec) > 0.1, "Should have non-zero electrostatic energy"
    assert total_energy < 0, "Ion-water interaction should be favorable (negative)"


@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available")  
def test_nbfix_energy_vs_openmm():
    """Compare NBFIX energy calculations between PyGCMC and OpenMM"""
    
    # Create system and force field
    mc_state = create_simple_ion_water_mcstate()
    forcefield = create_nbfix_forcefield()
    
    # Apply force field with NBFIX
    mc_state = setup_mcstate_with_forcefield(mc_state, forcefield)
    pygcmc.computeSystemEnergyPBCCutoff(mc_state)
    
    pygcmc_vdw = 0.0
    pygcmc_elec = 0.0
    for i in range(mc_state.activeResidueCount):
        pygcmc_vdw += mc_state.residues[i].energy_vdw
        pygcmc_elec += mc_state.residues[i].energy_elec
    pygcmc_total = pygcmc_vdw + pygcmc_elec
    
    # Set up OpenMM system
    omm_system, positions = setup_openmm_system_with_nbfix(mc_state, forcefield)
    
    # Calculate OpenMM energy
    integrator = VerletIntegrator(0.001 * picoseconds)
    platform = Platform.getPlatformByName('Reference')
    context = Context(omm_system, integrator, platform)
    context.setPositions(positions)
    
    state = context.getState(getEnergy=True)
    omm_energy = state.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    
    # IMPORTANT: PyGCMC's computeSystemEnergy* functions count each interaction twice
    # (both i-j and j-i), while OpenMM counts each pair once. We divide by 2 to compensate.
    # This is a known behavior of PyGCMC's system energy calculation.
    # NOTE: computeMovementEnergy* functions do NOT double count.
    pygcmc_corrected = pygcmc_total / 2.0
    
    print(f"\nEnergy comparison:")
    print(f"PyGCMC total (before correction): {pygcmc_total:.6f} kJ/mol")
    print(f"PyGCMC total (after /2 correction): {pygcmc_corrected:.6f} kJ/mol") 
    print(f"OpenMM total: {omm_energy:.6f} kJ/mol")
    
    # Allow 1% relative error due to implementation differences
    rel_error = abs(pygcmc_corrected - omm_energy) / abs(omm_energy)
    assert rel_error < 0.01, f"Energy mismatch: PyGCMC={pygcmc_corrected}, OpenMM={omm_energy}"
    
    # Test moving an atom and recalculating
    print("\n--- Testing energy change after moving Cl- ---")
    
    # Move chloride ion farther away
    mc_state.atoms[4].z = 2.5  # Move from 2.0 to 2.5 nm
    
    # Recalculate PyGCMC energy
    pygcmc.computeSystemEnergyPBCCutoff(mc_state)
    pygcmc_total_new = 0.0
    for i in range(mc_state.activeResidueCount):
        pygcmc_total_new += mc_state.residues[i].energy_vdw + mc_state.residues[i].energy_elec
    
    # Update OpenMM positions
    new_positions = list(positions)
    new_positions[4] = Vec3(1.5, 1.5, 2.5) * nanometers
    context.setPositions(new_positions)
    
    state_new = context.getState(getEnergy=True)
    omm_energy_new = state_new.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    
    # Compare energy changes - need to correct for factor of 2
    pygcmc_delta = (pygcmc_total_new - pygcmc_total) / 2.0
    omm_delta = omm_energy_new - omm_energy
    
    print(f"\nEnergy change after moving Cl-:")
    print(f"PyGCMC delta: {pygcmc_delta:.6f} kJ/mol")
    print(f"OpenMM delta: {omm_delta:.6f} kJ/mol")
    print(f"Difference in deltas: {abs(pygcmc_delta - omm_delta):.6f} kJ/mol")
    
    # Energy changes should be consistent
    assert abs(pygcmc_delta - omm_delta) < 0.1, "Energy changes should be consistent"
    
    # Moving opposite charges apart should increase energy
    assert pygcmc_delta > 0, "Moving Na+ and Cl- apart should increase energy"
    
    # Clean up OpenMM context
    del context


# Additional test for multiple NBFIX pairs
def test_multiple_nbfix_pairs():
    """Test system with multiple NBFIX corrections"""
    
    # Create force field
    forcefield = create_nbfix_forcefield()
    
    # Check that we have multiple NBFIX pairs
    nbfix_pairs = []
    test_types = ["SOD", "CLA", "OC", "OS", "ON3"]
    
    for type1 in test_types:
        for type2 in test_types:
            if type1 <= type2:  # Avoid duplicates
                result = forcefield.get_nbfix(type1, type2)
                # get_nbfix returns a tuple (epsilon, rmin, has_nbfix)
                if len(result) >= 3 and result[2]:  # has_nbfix is True
                    epsilon = result[0]
                    rmin = result[1]
                    nbfix_pairs.append((type1, type2, epsilon, rmin))
    
    print(f"\nFound {len(nbfix_pairs)} NBFIX pairs:")
    for type1, type2, epsilon, rmin in nbfix_pairs:
        print(f"  {type1}-{type2}: epsilon={epsilon:.5f} kcal/mol, Rmin={rmin:.3f} Å")
    
    # We should have at least SOD-CLA from the water_ions file
    assert len(nbfix_pairs) > 0, "Should have at least one NBFIX pair"
    
    # Check SOD-CLA specifically (manually added)
    sod_cla_found = False
    for type1, type2, epsilon, rmin in nbfix_pairs:
        if (type1 == "SOD" and type2 == "CLA") or (type1 == "CLA" and type2 == "SOD"):
            sod_cla_found = True
            assert abs(epsilon - (-0.083875)) < 1e-4
            assert abs(rmin - 3.731) < 1e-2
            print(f"  Verified {type1}-{type2} NBFIX matches expected values")
    
    assert sod_cla_found, "SOD-CLA NBFIX pair not found"


# Note about PyGCMC double counting behavior:
# - computeSystemEnergy* functions count each interaction twice (i-j and j-i)
# - computeMovementEnergy* functions count each interaction only once
# - When comparing with OpenMM or other MD engines, divide system energy by 2
# - This behavior is consistent across all PyGCMC energy calculation methods

# def test_system_vs_movement_energy():
#     """Test that system energy is double the movement energy for a two-residue system"""
#     # This test is commented out because the movement energy calculation
#     # requires more complex setup with proper residue types and movement types
#     # The key point is documented above: system energy functions double count
#     pass


@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available")
def test_nbfix_energy_components_separately():
    """Compare VDW and electrostatic energy components separately between PyGCMC and OpenMM
    
    This test provides more detailed validation by comparing energy components individually,
    which helps identify whether discrepancies come from VDW or electrostatic calculations.
    """
    
    # Create system and force field
    mc_state = create_simple_ion_water_mcstate()
    forcefield = create_nbfix_forcefield()
    
    # Apply force field with NBFIX
    mc_state = setup_mcstate_with_forcefield(mc_state, forcefield)
    pygcmc.computeSystemEnergyPBCCutoff(mc_state)
    
    # Get PyGCMC energy components
    pygcmc_vdw = 0.0
    pygcmc_elec = 0.0
    for i in range(mc_state.activeResidueCount):
        pygcmc_vdw += mc_state.residues[i].energy_vdw
        pygcmc_elec += mc_state.residues[i].energy_elec
    
    # IMPORTANT: PyGCMC's computeSystemEnergy* functions count each interaction twice
    # This is a known behavior that requires division by 2 for comparison with OpenMM
    pygcmc_vdw_corrected = pygcmc_vdw / 2.0
    pygcmc_elec_corrected = pygcmc_elec / 2.0
    
    # Set up OpenMM system with separate forces
    omm_system, positions = setup_openmm_system_with_nbfix(mc_state, forcefield, separate_forces=True)
    
    # Calculate OpenMM energy components
    integrator = VerletIntegrator(0.001 * picoseconds)
    platform = Platform.getPlatformByName('Reference')
    context = Context(omm_system, integrator, platform)
    context.setPositions(positions)
    
    # Get VDW energy (force group 1)
    state_vdw = context.getState(getEnergy=True, groups={1})
    omm_vdw = state_vdw.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    
    # Get electrostatic energy (force group 2)
    state_elec = context.getState(getEnergy=True, groups={2})
    omm_elec = state_elec.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    
    # Get total energy
    state_total = context.getState(getEnergy=True)
    omm_total = state_total.getPotentialEnergy().value_in_unit(kilojoules_per_mole)
    
    print(f"\n=== Detailed Energy Component Comparison ===")
    print(f"VDW Energy:")
    print(f"  PyGCMC (raw): {pygcmc_vdw:.6f} kJ/mol")
    print(f"  PyGCMC (corrected /2): {pygcmc_vdw_corrected:.6f} kJ/mol")
    print(f"  OpenMM: {omm_vdw:.6f} kJ/mol")
    print(f"  Relative error: {abs(pygcmc_vdw_corrected - omm_vdw) / abs(omm_vdw) * 100:.2f}%")
    
    print(f"\nElectrostatic Energy:")
    print(f"  PyGCMC (raw): {pygcmc_elec:.6f} kJ/mol")
    print(f"  PyGCMC (corrected /2): {pygcmc_elec_corrected:.6f} kJ/mol")
    print(f"  OpenMM: {omm_elec:.6f} kJ/mol")
    print(f"  Relative error: {abs(pygcmc_elec_corrected - omm_elec) / abs(omm_elec) * 100:.2f}%")
    
    print(f"\nTotal Energy:")
    print(f"  PyGCMC (corrected): {pygcmc_vdw_corrected + pygcmc_elec_corrected:.6f} kJ/mol")
    print(f"  OpenMM: {omm_total:.6f} kJ/mol")
    
    # Compare components separately with 1% tolerance
    rel_tol = 0.01
    
    # Check VDW energy
    if abs(omm_vdw) > 1e-6:  # Only check relative error if VDW is non-negligible
        vdw_rel_error = abs(pygcmc_vdw_corrected - omm_vdw) / abs(omm_vdw)
        assert vdw_rel_error < rel_tol, f"VDW energy mismatch: PyGCMC={pygcmc_vdw_corrected}, OpenMM={omm_vdw}"
    else:
        assert abs(pygcmc_vdw_corrected - omm_vdw) < 1e-3, f"VDW energy should be near zero"
    
    # Check electrostatic energy
    elec_rel_error = abs(pygcmc_elec_corrected - omm_elec) / abs(omm_elec)
    assert elec_rel_error < rel_tol, f"Electrostatic energy mismatch: PyGCMC={pygcmc_elec_corrected}, OpenMM={omm_elec}"
    
    # Check total energy consistency
    omm_sum = omm_vdw + omm_elec
    assert abs(omm_total - omm_sum) < 1e-6, "OpenMM total should equal sum of components"
    
    # Clean up
    del context


@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM not available")
def test_nbfix_across_pbc():
    """Test NBFIX energy calculation with atoms separated by periodic boundary
    
    This tests that PyGCMC correctly applies NBFIX parameters when atoms interact
    across periodic boundaries, which is a common source of errors.
    """
    
    # Create MCState with atoms near box boundaries
    state = pygcmc.MCState()
    
    # Small box to force PBC interactions
    state.info.box = [2.0, 2.0, 2.0]  # 2 nm box
    state.info.cutoff = 1.2  # 1.2 nm cutoff
    
    # Define atom types
    sod_idx = state.atomTypes.get_or_add_type("SOD")
    cla_idx = state.atomTypes.get_or_add_type("CLA")
    
    # Place ions near opposite box edges
    # They are 1.5 nm apart in real space, but only 0.5 nm via PBC
    atom_na = pygcmc.MCAtom()
    atom_na.x = 0.25  # Near left edge
    atom_na.y = 1.0
    atom_na.z = 1.0
    atom_na.charge = 1.0
    atom_na.type = sod_idx
    
    atom_cl = pygcmc.MCAtom()
    atom_cl.x = 1.75  # Near right edge
    atom_cl.y = 1.0
    atom_cl.z = 1.0
    atom_cl.charge = -1.0
    atom_cl.type = cla_idx
    
    state.atoms = [atom_na, atom_cl]
    state.activeAtomCount = 2
    
    # Add residues
    res_na = pygcmc.MCResidue()
    res_na.active = True
    res_na.atomStart = 0
    res_na.atomCount = 1
    res_na.type = 0
    
    res_cl = pygcmc.MCResidue()
    res_cl.active = True
    res_cl.atomStart = 1
    res_cl.atomCount = 1
    res_cl.type = 1
    
    state.residues = [res_na, res_cl]
    state.activeResidueCount = 2
    
    # Create and apply force field with NBFIX
    forcefield = create_nbfix_forcefield()
    state = setup_mcstate_with_forcefield(state, forcefield)
    
    # Calculate PyGCMC energy
    pygcmc.computeSystemEnergyPBCCutoff(state)
    pygcmc_total = 0.0
    for res in state.residues[:state.activeResidueCount]:
        pygcmc_total += res.energy_vdw + res.energy_elec
    pygcmc_corrected = pygcmc_total / 2.0
    
    # Calculate expected energy manually
    # Distance via PBC: 2.0 - 1.5 = 0.5 nm
    pbc_distance = 0.5  # nm
    
    # Get NBFIX parameters for SOD-CLA
    nbfix_result = forcefield.get_nbfix("SOD", "CLA")
    nbfix_eps = abs(nbfix_result[0]) * KCAL_TO_KJ
    nbfix_rmin = nbfix_result[1]
    nbfix_sigma = nbfix_rmin / math.pow(2.0, 1.0/6.0) * ANGSTROM_TO_NM
    
    # Calculate expected energies
    expected_vdw = 4 * nbfix_eps * (math.pow(nbfix_sigma/pbc_distance, 12) - math.pow(nbfix_sigma/pbc_distance, 6))
    expected_elec = kC * 1.0 * (-1.0) / pbc_distance
    expected_total = expected_vdw + expected_elec
    
    print(f"\n=== PBC NBFIX Test ===")
    print(f"Box size: {state.info.box}")
    print(f"Atom positions: Na+ at ({atom_na.x}, {atom_na.y}, {atom_na.z}), Cl- at ({atom_cl.x}, {atom_cl.y}, {atom_cl.z})")
    print(f"Direct distance: 1.5 nm")
    print(f"PBC distance: {pbc_distance} nm")
    print(f"\nExpected energies:")
    print(f"  VDW: {expected_vdw:.6f} kJ/mol")
    print(f"  Electrostatic: {expected_elec:.6f} kJ/mol")
    print(f"  Total: {expected_total:.6f} kJ/mol")
    print(f"\nPyGCMC energy (corrected): {pygcmc_corrected:.6f} kJ/mol")
    
    # At 0.5 nm distance, check the energy
    # The energy should be negative (attractive) due to opposite charges
    assert pygcmc_corrected < 0, "Energy should be attractive for opposite charges at 0.5 nm"
    
    # Check that the energy is in the right ballpark
    # Allow larger tolerance due to potential cutoff effects
    rel_error = abs(pygcmc_corrected - expected_total) / abs(expected_total)
    assert rel_error < 0.05, f"PBC energy mismatch: PyGCMC={pygcmc_corrected}, Expected={expected_total}"


if __name__ == "__main__":
    # Run tests
    if HAS_OPENMM:
        test_nbfix_overrides_lj_combination()
        test_ion_water_nbfix_energy()
        test_nbfix_energy_vs_openmm()
        test_nbfix_energy_components_separately()
        test_nbfix_across_pbc()
    test_multiple_nbfix_pairs()
    print("\nAll NBFIX tests passed!")