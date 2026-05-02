#pragma once

// Include all energy calculation modules
#include "common/EnergyMain.hpp"
#include "lj/LJMain.hpp"
#include "coulomb/CoulombMain.hpp"
#include "ewald/EwaldMain.hpp"
#include "pme/PMEMain.hpp"

/**
 * @file EnergyModule.hpp
 * @brief Energy Module - Unified Entry Point for All Energy Calculation Functionality
 *
 * This is the ONLY header file you need to include to access all energy calculation
 * capabilities in the GCMC simulation framework. The module provides multiple energy
 * calculation methods optimized for different system sizes and accuracy requirements.
 *
 * **Module Organization:**
 *
 * **common/ directory** - Core Energy Infrastructure
 * - EnergyMain.hpp: Main interface with unified energy calculation functions
 * - EnergyInterface.hpp: Abstract interfaces and enums (EnergyMethod, unified API)
 * - EnergyDirectCore.hpp: Direct summation for nonbonded interactions
 * - EnergyUtils.hpp: Utility functions for energy calculations
 * - EnergyConstants.hpp: Physical constants and conversion factors
 *
 * **lj/ directory** - Lennard-Jones (van der Waals) Interactions
 * - LJMain.hpp: Main interface for LJ energy calculations
 * - LJPotential.hpp: Core Lennard-Jones potential implementations
 * - LJSwitch.hpp: CHARMM-style smooth switching functions for LJ
 * - LJForceField.hpp: Force field parameter management for LJ interactions
 *
 * **coulomb/ directory** - Coulombic (Electrostatic) Interactions
 * - CoulombMain.hpp: Main interface for Coulombic energy calculations
 * - CoulombPairCore.hpp: Pair-wise Coulombic interaction calculations
 * - CoulombPotential.hpp: Core Coulombic potential implementations
 * - CoulombUtils.hpp: Utility functions for electrostatic calculations
 *
 * **ewald/ directory** - Ewald Summation Method
 * - EwaldMain.hpp: Main interface for Ewald summation calculations
 * - EwaldComposite.hpp: High-level Ewald calculation coordination
 * - EwaldCore.hpp: Core Ewald summation algorithm implementation
 * - EwaldReal.hpp: Real-space contribution calculations
 * - EwaldRecip.hpp: Reciprocal-space contribution calculations (k-space)
 * - EwaldSelf.hpp: Self-energy correction calculations
 * - EwaldUtils.hpp: Utility functions for Ewald calculations
 *
 * **pme/ directory** - Particle Mesh Ewald Method
 * - PMEMain.hpp: Main interface for PME calculations
 * - PMEComposite.hpp: High-level PME calculation coordination
 * - PMECore.hpp: Core PME algorithm implementation
 * - PMEGrid.hpp: Grid-based charge assignment and force interpolation
 * - PMEFFTCore.hpp: Core FFT algorithms
 * - PMEFFT3D.hpp: 3D FFT batch operations for PME
 * - PMEBSpline.hpp: B-spline interpolation for charge assignment
 * - PMEReal.hpp: Real-space contribution calculations (similar to Ewald)
 * - PMERecip.hpp: Reciprocal-space calculations using FFT
 * - PMESelf.hpp: Self-energy corrections for PME
 * - PMEUtils.hpp: Utility functions and parameter optimization
 *
 * **Energy Calculation Methods:**
 *
 * 1. **DIRECT**: Explicit pair-wise summation
 *    - Fast for small systems (< 1000 atoms)
 *    - O(N²) scaling
 *    - Supports cutoff and PBC options
 *
 * 2. **EWALD**: Ewald summation for long-range electrostatics
 *    - Optimal for medium systems (1000-10000 atoms)
 *    - O(N^1.5) scaling
 *    - Exact treatment of long-range interactions in periodic systems
 *
 * 3. **PME**: Particle Mesh Ewald with FFT acceleration
 *    - Best for large systems (> 10000 atoms)
 *    - O(N log N) scaling
 *    - Fast approximation of Ewald summation using grid-based FFT
 *
 * **Usage Examples:**
 *
 * ```cpp
 * #include "platform/cpu/energy/EnergyModule.hpp"
 * using namespace pygcmc::platform::cpu;
 *
 * // 1. Direct calculation for small systems
 * computeSystemEnergyDirect(state, false, false);  // No cutoff, no PBC
 * computeSystemEnergyDirect(state, true, true);    // With cutoff and PBC
 *
 * // 2. Using unified interface with method selection
 * computeSystemEnergy(state, EnergyMethod::DIRECT, true, true);
 * computeSystemEnergy(state, EnergyMethod::EWALD);
 * computeSystemEnergy(state, EnergyMethod::PME);
 *
 * // 3. Initialize Ewald parameters before use
 * initializeEwaldParameters(cutoff, box);
 * double ewald_energy = energy::getTotalEnergy(state, EnergyMethod::EWALD);
 *
 * // 4. Initialize PME parameters before use
 * initializePMEParameters(cutoff, box, grid_size);
 * double pme_energy = energy::getTotalEnergy(state, EnergyMethod::PME);
 *
 * // 5. Get method information
 * std::string method_name = energy::getEnergyMethodName(EnergyMethod::PME);
 * std::cout << "Using method: " << method_name << std::endl;
 * ```
 *
 * **Performance Guidelines:**
 *
 * - **< 1,000 atoms**: Use DIRECT method for simplicity and speed
 * - **1,000 - 10,000 atoms**: Use EWALD for accuracy with reasonable performance
 * - **> 10,000 atoms**: Use PME for optimal performance in large systems
 * - **Cutoff radius**: Typically 1.0-1.4 nm for good accuracy/performance balance
 * - **Grid spacing**: 0.1-0.12 nm for PME calculations provides good accuracy
 *
 * **Design Philosophy:**
 * - Single header inclusion for all energy functionality
 * - Method-agnostic unified interface for easy switching
 * - Optimized implementations for different system sizes
 * - Well-structured: each file under 300 lines with focused responsibilities
 * - Performance-critical: minimal overhead through careful design
 * - Backward compatible: maintains existing API while providing new features
 */
