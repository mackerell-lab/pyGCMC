#pragma once

#ifndef PYGCMC_MODEL_MODULE_HPP
#define PYGCMC_MODEL_MODULE_HPP

/**
 * @file ModelModule.hpp
 * @brief Model Module - Data Layer Unified Entry Point
 * 
 * This is the single entry point for the entire Model module, providing access to
 * all data structures and utilities needed for molecular modeling and GCMC simulation.
 * 
 * The Model module is organized into several functional areas:
 * 
 * 📁 common/     - Common interfaces, constants, and utilities
 * 📁 atom/       - Atom-level data structures and operations  
 * 📁 residue/    - Residue-level composition and validation
 * 📁 molecule/   - Molecular system management and utilities
 * 📁 topology/   - Topology and force field parameters
 * 📁 montecarlo/ - Monte Carlo state and simulation data
 * 📁 param/      - Global simulation parameters
 * 
 * Design Principles:
 * - AI-Friendly: Small, focused files (≤200 lines each)
 * - Maximum Extensibility: Clear separation of concerns
 * - Consistent Naming: Main/Composite/Core/Utils pattern
 * - Backward Compatibility: All original APIs preserved
 */

// === Core Interfaces ===
#include "common/ModelInterface.hpp"
#include "common/ModelConstants.hpp"
#include "common/ModelUtils.hpp"

// === Atom Layer ===
#include "atom/AtomMain.hpp"

// === Residue Layer ===
#include "residue/ResidueMain.hpp"

// === Molecule Layer ===
#include "molecule/MolecularMain.hpp"

// === Topology & Force Field ===
#include "topology/TopologyMain.hpp"
#include "topology/ForceFieldMain.hpp"

// === Monte Carlo Simulation ===
#include "montecarlo/MCMain.hpp"

// === Parameters ===
#include "param/ParamMain.hpp"

// === Legacy Structure Support ===
#include "structure.hpp"

namespace pygcmc {
namespace model {

/**
 * @brief Model Module Information
 */
namespace info {
    constexpr const char* VERSION = "2.0.0";
    constexpr const char* BUILD_DATE = __DATE__;
    
    constexpr int TOTAL_COMPONENTS = 7;
    constexpr const char* COMPONENTS[] = {
        "common", "atom", "residue", "molecule", 
        "topology", "montecarlo", "param"
    };
}

/**
 * @brief Quick Start Example: Creating a Simple Molecular System
 * 
 * @code{.cpp}
 * #include "model/ModelModule.hpp"
 * using namespace pygcmc::model;
 * 
 * // 1. Create atoms
 * auto atom1 = std::make_shared<Atom>("CA", "CT1", 0.07, 12.01, 
 *                                     1.0, 2.0, 3.0, "ALA", 1, "PROT");
 * auto atom2 = std::make_shared<Atom>("CB", "CT3", -0.27, 12.01,
 *                                     2.0, 3.0, 4.0, "ALA", 1, "PROT");
 * 
 * // 2. Create residue and add atoms
 * auto residue = std::make_shared<Residue>("ALA", 1, "PROT");
 * residue->add_atom(atom1);
 * residue->add_atom(atom2);
 * residue->calculate_center_of_mass();
 * 
 * // 3. Create molecular system
 * auto molecular = std::make_shared<Molecular>();
 * molecular->add_atom(atom1);
 * molecular->add_atom(atom2);
 * molecular->add_residue(residue);
 * 
 * // 4. Set up topology
 * auto topology = std::make_shared<Topology>();
 * int atom1_id = topology->add_atom("CA", "CT1", 0.07, 12.01, "ALA", 1, "PROT");
 * int atom2_id = topology->add_atom("CB", "CT3", -0.27, 12.01, "ALA", 1, "PROT");
 * topology->add_bond(atom1_id, atom2_id, 1.53, 340.0);
 * 
 * // 5. Set up force field
 * auto forcefield = std::make_shared<ForceField>();
 * forcefield->add_atom_mass("CT1", 12.01);
 * forcefield->add_atom_mass("CT3", 12.01);
 * forcefield->add_lj_params("CT1", -0.020, 2.275);
 * forcefield->add_lj_params("CT3", -0.080, 2.060);
 * 
 * // 6. Set up MC state (for GCMC simulation)
 * auto mc_state = std::make_shared<MCState>();
 * mc_state->setTemperature(300.0);
 * mc_state->setBoxDimensions(50.0, 50.0, 50.0);
 * 
 * // 7. Set up parameters
 * auto params = std::make_shared<Param>();
 * params->set_temperature(300.0);
 * params->set_box_size(50.0, 50.0, 50.0);
 * params->set_mc_steps(1000000);
 * @endcode
 */

/**
 * @brief Advanced Example: GCMC Simulation Setup
 * 
 * @code{.cpp}
 * #include "model/ModelModule.hpp"
 * using namespace pygcmc::model;
 * 
 * // Create a complete GCMC system
 * void setup_gcmc_simulation() {
 *     // 1. Initialize parameters
 *     auto params = std::make_shared<Param>();
 *     auto& mc_info = params->get_mc_info();
 *     mc_info.setTemperature(300.0);
 *     mc_info.mc_steps = 1000000;
 *     mc_info.insertion_deletion_frac = 0.5;
 * 
 *     auto& space_info = params->get_space_info();
 *     space_info.setBoxSize(50.0, 50.0, 50.0);
 *     space_info.cutoff = 12.0;
 * 
 *     auto& fragment_info = params->get_fragment_info();
 *     fragment_info.water_density = 55.0;  // M
 *     fragment_info.conc_list = {0.25, 55.0};  // solute and water concentrations
 *     fragment_info.muex_list = {-5.0, -6.1};  // excess chemical potentials
 * 
 *     // 2. Set up MC state
 *     auto mc_state = std::make_shared<MCState>();
 *     mc_state->setTemperature(300.0);
 *     mc_state->setBoxDimensions(50.0, 50.0, 50.0);
 *     mc_state->setCutoff(12.0);
 * 
 *     // 3. Set up force field for water and solute
 *     mc_state->setupForceField(4);  // 4 atom types: O, H, C, N
 *     mc_state->setLJParameters(0, 0, 3.15, 0.636);  // O-O
 *     mc_state->setLJParameters(1, 1, 1.00, 0.046);  // H-H
 *     mc_state->setLJParameters(0, 1, 2.08, 0.192);  // O-H (mixed)
 * 
 *     // 4. Add water molecules
 *     int water_type = mc_state->getOrAddResidueType("SOL");
 *     mc_state->addMovementResidue("SOL", 0, 1000);  // up to 1000 water molecules
 * 
 *     // System is now ready for GCMC simulation
 * }
 * @endcode
 */

/**
 * @brief Utility Functions for Common Operations
 */
namespace utils {

/**
 * @brief Create a water molecule
 */
inline std::shared_ptr<Residue> create_water_molecule(int residue_id, const std::string& segment = "SOLV") {
    auto water = std::make_shared<Residue>("SOL", residue_id, segment);
    
    // Add atoms: O, H1, H2
    auto oxygen = std::make_shared<Atom>("O", "OT", -0.834, 15.999, 
                                         0.0, 0.0, 0.0, "SOL", residue_id, segment);
    auto h1 = std::make_shared<Atom>("H1", "HT", 0.417, 1.008,
                                     0.757, 0.586, 0.0, "SOL", residue_id, segment);
    auto h2 = std::make_shared<Atom>("H2", "HT", 0.417, 1.008,
                                     -0.757, 0.586, 0.0, "SOL", residue_id, segment);
    
    water->add_atom(oxygen);
    water->add_atom(h1);
    water->add_atom(h2);
    water->calculate_center_of_mass();
    
    return water;
}

/**
 * @brief Create a simple protein residue (Alanine)
 */
inline std::shared_ptr<Residue> create_alanine_residue(int residue_id, const std::string& segment = "PROT") {
    auto ala = std::make_shared<Residue>("ALA", residue_id, segment);
    
    // Add backbone atoms
    auto n = std::make_shared<Atom>("N", "NH1", -0.47, 14.007,
                                    0.0, 0.0, 0.0, "ALA", residue_id, segment);
    auto ca = std::make_shared<Atom>("CA", "CT1", 0.07, 12.01,
                                     1.458, 0.0, 0.0, "ALA", residue_id, segment);
    auto c = std::make_shared<Atom>("C", "C", 0.51, 12.01,
                                    2.009, 1.421, 0.0, "ALA", residue_id, segment);
    auto o = std::make_shared<Atom>("O", "O", -0.51, 15.999,
                                    1.239, 2.364, 0.0, "ALA", residue_id, segment);
    auto cb = std::make_shared<Atom>("CB", "CT3", -0.27, 12.01,
                                     2.196, -0.889, -1.07, "ALA", residue_id, segment);
    
    ala->add_atom(n);
    ala->add_atom(ca);
    ala->add_atom(c);
    ala->add_atom(o);
    ala->add_atom(cb);
    ala->calculate_center_of_mass();
    
    return ala;
}

/**
 * @brief Validate a complete molecular system
 */
inline bool validate_molecular_system(const Molecular& system) {
    // Check basic validity
    if (!system.is_valid()) return false;
    
    // Check all residues
    for (const auto& residue : system.get_residues()) {
        if (!residue || !residue->is_valid()) return false;
    }
    
    // Check all atoms
    for (const auto& atom : system.get_atoms()) {
        if (!atom || !atom->is_valid()) return false;
    }
    
    return true;
}

/**
 * @brief Get system statistics
 */
inline std::string get_system_summary(const Molecular& system) {
    auto stats = system.get_system_statistics();
    std::stringstream ss;
    ss << "System Summary:\n";
    ss << "  Atoms: " << stats.num_atoms << " (" << stats.num_heavy_atoms 
       << " heavy, " << stats.num_hydrogen_atoms << " H)\n";
    ss << "  Residues: " << stats.num_residues << "\n";
    ss << "  Chains: " << stats.num_chains << "\n";
    ss << "  Segments: " << stats.num_segments << "\n";
    ss << "  Total Mass: " << stats.total_mass << " amu\n";
    ss << "  Total Charge: " << stats.total_charge << " e\n";
    ss << "  Center of Mass: [" << stats.center_of_mass[0] << ", " 
       << stats.center_of_mass[1] << ", " << stats.center_of_mass[2] << "]\n";
    return ss.str();
}

} // namespace utils

/**
 * @brief Module Testing and Validation
 */
namespace test {

/**
 * @brief Run basic functionality tests for all components
 */
inline bool run_basic_tests() {
    try {
        // Test atom creation
        auto atom = std::make_shared<Atom>("CA", "CT1", 0.07, 12.01, 0.0, 0.0, 0.0, "ALA", 1, "PROT");
        if (!atom->is_valid()) return false;
        
        // Test residue creation
        auto residue = utils::create_alanine_residue(1);
        if (!residue->is_valid()) return false;
        
        // Test molecular system
        auto molecular = std::make_shared<Molecular>();
        molecular->add_residue(residue);
        if (!molecular->is_valid()) return false;
        
        // Test topology
        auto topology = std::make_shared<Topology>();
        topology->add_atom("CA", "CT1", 0.07, 12.01, "ALA", 1, "PROT");
        if (!topology->is_valid()) return false;
        
        // Test force field
        auto ff = std::make_shared<ForceField>();
        ff->add_atom_mass("CT1", 12.01);
        ff->add_lj_params("CT1", 0.02, 2.275);
        if (!ff->is_valid()) return false;
        
        // Test MC state
        auto mc = std::make_shared<MCState>();
        mc->setTemperature(300.0);
        if (!mc->is_valid()) return false;
        
        // Test parameters
        auto params = std::make_shared<Param>();
        params->set_temperature(300.0);
        if (!params->is_valid()) return false;
        
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

} // namespace test

} // namespace model
} // namespace pygcmc

/**
 * @brief Module Version and Build Information
 */
#define PYGCMC_MODEL_VERSION_MAJOR 2
#define PYGCMC_MODEL_VERSION_MINOR 0
#define PYGCMC_MODEL_VERSION_PATCH 0
#define PYGCMC_MODEL_VERSION_STRING "2.0.0"

#endif // PYGCMC_MODEL_MODULE_HPP 