# GCMC Movement Module Tests

This directory contains comprehensive tests for the GCMC (Grand Canonical Monte Carlo) movement module.

## Test Structure

The tests are organized into modular files, each under 300 lines:

### Main Test File
- `test_movement_gcmc.py` - Main entry point that imports all sub-tests

### Test Modules

1. **basic_operations.py** (89 lines)
   - `test_insertion_attempt` - Test basic molecule insertion
   - `test_deletion_attempt` - Test molecule deletion
   - `test_translation_attempt` - Test molecule translation
   - `test_rotation_attempt` - Test molecule rotation

2. **acceptance_tests.py** (108 lines)
   - `test_acceptance_rate_range` - Verify acceptance rates are in [0,1]
   - `test_chemical_potential_effect` - Test μ effect on insertions
   - `test_temperature_effect` - Test temperature effect on acceptance

3. **energy_tests.py** (69 lines)
   - `test_energy_conservation` - Verify energy conservation in rejected moves
   - `test_fragment_energy_calculation` - Test fragment energy calculations
   - `test_system_energy_consistency` - Check energy calculation consistency

4. **cavity_bias_tests.py** (125 lines)
   - `test_cavity_detection` - Test cavity finding in empty box
   - `test_cavity_bias_insertion` - Test biased vs uniform insertion
   - `test_cavity_grid_parameters` - Test grid spacing effects
   - `test_cavity_vs_uniform_insertion` - Compare insertion methods

5. **reservoir_tests.py** (152 lines)
   - `test_template_management` - Test fragment template operations
   - `test_instance_creation_deletion` - Test instance lifecycle
   - `test_ghost_fragment_recycling` - Test ghost recycling mechanism
   - `test_reservoir_statistics` - Test statistics tracking

6. **scaling_tests.py** (186 lines)
   - `test_volume_scaling` - Test N ∝ V relationship
   - `test_chemical_potential_scaling` - Test N ∝ exp(βμ) relationship
   - `test_detailed_balance` - Verify detailed balance in equilibrium

### Fixtures (conftest.py)

Common fixtures used across all tests:
- `setup_gcmc_state` - Create basic MCState with water parameters
- `gcmc_params` - Standard GCMC movement parameters
- `water_template` - TIP3P water molecule template
- `gcmc_mover` - Configured MovementModule
- `small_state` - Small box for testing
- `populated_state` - State with pre-inserted molecules

Helper functions:
- `calculate_acceptance_rate()` - Calculate overall acceptance
- `run_gcmc_steps()` - Run specified number of GCMC steps

## Running Tests

Run all GCMC tests:
```bash
pytest tests/simulation/test_movement_gcmc.py
```

Run specific test module:
```bash
pytest tests/simulation/movementGCMC/basic_operations.py
```

Run single test:
```bash
pytest tests/simulation/movementGCMC/basic_operations.py::test_insertion_attempt -v
```

## Key Concepts Tested

1. **Grand Canonical Ensemble (μVT)**
   - Fixed chemical potential, volume, temperature
   - Variable particle number

2. **Move Types**
   - Insertion: Add molecule from reservoir
   - Deletion: Remove molecule to reservoir  
   - Translation: Displace molecule
   - Rotation: Rotate molecule orientation

3. **Acceptance Criteria**
   - Metropolis criterion with proper μVT factors
   - Detailed balance verification

4. **Biasing Techniques**
   - Cavity bias: Insert in favorable locations
   - Configurational bias: Multi-trial insertions

5. **Energy Calculations**
   - Van der Waals (Lennard-Jones)
   - Electrostatic (if charges present)
   - Proper cutoffs and PBC

## Dependencies

- PyGCMC C++ bindings
- NumPy for numerical operations
- pytest for testing framework

## Notes

- Some tests may be skipped if certain features (like FragmentReservoir) are not exposed in Python bindings
- Test parameters are tuned for TIP3P water at 300K
- Acceptance rates and scaling tests have wide tolerances due to Monte Carlo noise