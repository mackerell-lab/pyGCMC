#pragma once

// Include all system main modules
#include "common/SystemMain.hpp"
#include "common/SystemFactory.hpp"
#include "log/LogMain.hpp"
#include "molecular/MolecularMain.hpp"
#include "montecarlo/MCMain.hpp"

/**
 * @file SystemModule.hpp
 * @brief System Module - Unified Entry Point for All System Functionality
 *
 * This is the ONLY header file you need to include to access all system management
 * capabilities in the GCMC simulation framework. The module is organized into
 * specialized sub-modules, each handling specific aspects of system management.
 *
 * **Module Organization:**
 *
 * **common/ directory** - Core System Infrastructure
 * - SystemMain.hpp: Original System class for parameter management and initialization
 * - SystemFactory.hpp: Factory functions and version information
 * - SystemInterface.hpp: Abstract interfaces and enums (LogLevel, SystemKind)
 * - SystemConstants.hpp: Physical constants and unit conversions
 * - SystemUtils.hpp: Utility functions (PBC, string manipulation, geometric calculations)
 *
 * **log/ directory** - Logging System
 * - LogMain.hpp: Comprehensive logging with configurable levels (DEBUG, INFO, WARNING, ERROR)
 * - Static interface: LogMain::set_verbose(), LogMain::set_log_level()
 * - Template logging: LogMain::log(), LogMain::debug(), LogMain::info()
 * - Convenience functions: initializeLogging(), setVerbose(), setLogLevel()
 *
 * **molecular/ directory** - Molecular System Management
 * - MolecularMain.hpp: Main interface (backward compatible MolecularSystem)
 * - MolecularComposite.hpp: High-level composition and coordination
 * - MolecularCombiner.hpp: Combine PDB Structure with Topology files
 * - MolecularMatcher.hpp: Match residue sequences between PDB and topology
 * - MolecularMerger.hpp: Merge multiple topology files into molecular system
 * - MolecularValidator.hpp: Validate atom types and topology consistency
 *
 * **montecarlo/ directory** - Monte Carlo GCMC Simulation Engine
 * - MCMain.hpp: Main interface (backward compatible MonteCarloSystem)
 * - MCComposite.hpp: High-level GCMC operations coordination
 * - MCCore.hpp: Core GCMC operations (insert, remove, translate residues)
 * - MCInitializer.hpp: Initialize GCMC system from molecular structure
 * - MCGeometry.hpp: Periodic boundary conditions and geometric calculations
 * - MCSwitching.hpp: CHARMM-style smooth switching functions
 * - MCMovementBuilder.hpp: Build movement molecule configurations
 * - MCMovementTypeCollector.hpp: Collect and manage movement molecule types
 * - MCMovementReindexer.hpp: Reindex existing residues for GCMC operations
 *
 * **Usage Examples:**
 *
 * ```cpp
 * #include "system/SystemModule.hpp"
 * using namespace pygcmc::system;
 *
 * // 1. Initialize logging
 * initializeLogging(true, LogLevel::DEBUG);
 *
 * // 2. Build molecular system from files
 * MolecularSystem molSys;
 * auto structure = io::readPDB("protein.pdb");
 * auto topology = io::readTopology("charmm36.top");
 * auto molecular = molSys.combine(structure, topology);
 *
 * // 3. Set up Monte Carlo simulation
 * MonteCarloSystem mcSys;
 * mcSys.initializeFromMolecular(molecular);
 * mcSys.setSwitchingFunction(true, 1.0f, 1.2f);
 *
 * // 4. Add movement molecules for GCMC
 * std::vector<MovementMolecularInfo> molecules = {
 *     {waterMolecular, 1000},  // Max 1000 water molecules
 *     {ionMolecular, 50}       // Max 50 ion pairs
 * };
 * mcSys.addMovementMolecules(molecules);
 *
 * // 5. Perform GCMC operations
 * int inserted = mcSys.insertResidue(waterResidue, waterAtoms);
 * float energy = mcSys.calcTotalEnergy();
 * bool removed = mcSys.removeResidue(inserted);
 * ```
 *
 * **Backward Compatibility:**
 * All original APIs (System, MolecularSystem, MonteCarloSystem) are preserved
 * through using declarations in each sub-module. Existing code requires no changes.
 *
 * **Design Philosophy:**
 * - Single header inclusion for all functionality
 * - Modular architecture with clear separation of concerns
 * - Each sub-module is independently testable and maintainable
 * - Well-structured: each file under 300 lines with focused responsibilities
 * - Performance-oriented: minimal overhead through careful design
 */
