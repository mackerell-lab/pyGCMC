#pragma once

#ifndef PYGCMC_MODEL_UTILS_FACTORY_HPP
#define PYGCMC_MODEL_UTILS_FACTORY_HPP

#include "../atom/AtomMain.hpp"
#include "../residue/ResidueMain.hpp"
#include "../molecule/MolecularMain.hpp"
#include <memory>
#include <string>
#include <sstream>

namespace pygcmc {
namespace model {
namespace factory {

// Using declarations for cleaner code
using AtomPtr = std::shared_ptr<atom::Atom>;
using ResiduePtr = std::shared_ptr<residue::Residue>;
using MolecularPtr = std::shared_ptr<molecule::Molecular>;

/**
 * @brief Create a water molecule
 */
inline ResiduePtr create_water_molecule(int residue_id, const std::string& segment = "SOLV") {
    auto water = std::make_shared<residue::Residue>("SOL", residue_id, segment);
    
    // Add atoms: O, H1, H2 with correct constructor parameters
    // AtomCore(int bynu, type, resname, ires, segid, iseg, x, y, z, wmain, mass, charge, chem, hetatm)
    auto oxygen = std::make_shared<atom::Atom>(1, "OT", "SOL", residue_id, segment, 0,
                                               0.0, 0.0, 0.0, 1.0, 15.999, -0.834, "O", false);
    auto h1 = std::make_shared<atom::Atom>(2, "HT", "SOL", residue_id, segment, 0,
                                           0.757, 0.586, 0.0, 1.0, 1.008, 0.417, "H", false);
    auto h2 = std::make_shared<atom::Atom>(3, "HT", "SOL", residue_id, segment, 0,
                                           -0.757, 0.586, 0.0, 1.0, 1.008, 0.417, "H", false);
    
    // Add atoms to residue
    water->add_atom(oxygen);
    water->add_atom(h1);
    water->add_atom(h2);
    
    return water;
}

/**
 * @brief Create a simple protein residue (Alanine)
 */
inline ResiduePtr create_alanine_residue(int residue_id, const std::string& segment = "PROT") {
    auto ala = std::make_shared<residue::Residue>("ALA", residue_id, segment);
    
    // Add backbone atoms with correct constructor parameters
    auto n = std::make_shared<atom::Atom>(1, "NH1", "ALA", residue_id, segment, 0,
                                          0.0, 0.0, 0.0, 1.0, 14.007, -0.47, "N", false);
    auto ca = std::make_shared<atom::Atom>(2, "CT1", "ALA", residue_id, segment, 0,
                                           1.458, 0.0, 0.0, 1.0, 12.01, 0.07, "C", false);
    auto c = std::make_shared<atom::Atom>(3, "C", "ALA", residue_id, segment, 0,
                                          2.009, 1.421, 0.0, 1.0, 12.01, 0.51, "C", false);
    auto o = std::make_shared<atom::Atom>(4, "O", "ALA", residue_id, segment, 0,
                                          1.239, 2.364, 0.0, 1.0, 15.999, -0.51, "O", false);
    auto cb = std::make_shared<atom::Atom>(5, "CT3", "ALA", residue_id, segment, 0,
                                           2.196, -0.889, -1.07, 1.0, 12.01, -0.27, "C", false);
    
    // Add atoms to residue
    ala->add_atom(n);
    ala->add_atom(ca);
    ala->add_atom(c);
    ala->add_atom(o);
    ala->add_atom(cb);
    
    return ala;
}

/**
 * @brief Validate a complete molecular system
 */
inline bool validate_molecular_system(const molecule::Molecular& system) {
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
inline std::string get_system_summary(const molecule::Molecular& system) {
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

/**
 * @brief Run basic functionality tests for all components
 */
inline bool run_basic_tests() {
    try {
        // Test atom creation
        auto atom = std::make_shared<atom::Atom>(1, "CT1", "ALA", 1, "PROT", 0, 
                                                  0.0, 0.0, 0.0, 1.0, 12.01, 0.07, "C", false);
        if (!atom->is_valid()) return false;
        
        // Test residue creation
        auto residue = create_alanine_residue(1);
        if (!residue->is_valid()) return false;
        
        // Test molecular system
        auto molecular = std::make_shared<molecule::Molecular>();
        molecular->add_residue(residue);
        if (!molecular->is_valid()) return false;
        
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

/**
 * @brief Test factory functions
 */
inline bool test_factory_functions() {
    try {
        // Test water molecule creation
        auto water = create_water_molecule(1);
        if (!water || !water->is_valid()) return false;
        
        // Test alanine residue creation
        auto ala = create_alanine_residue(2);
        if (!ala || !ala->is_valid()) return false;
        
        // Test molecular system validation
        auto molecular = std::make_shared<molecule::Molecular>();
        molecular->add_residue(water);
        molecular->add_residue(ala);
        
        if (!validate_molecular_system(*molecular)) return false;
        
        // Test system summary
        std::string summary = get_system_summary(*molecular);
        if (summary.empty()) return false;
        
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

/**
 * @brief Run all tests
 */
inline bool run_all_tests() {
    return run_basic_tests() && test_factory_functions();
}

} // namespace factory

/**
 * @brief Module Information and Version (following system pattern)
 */
namespace info {
    constexpr const char* VERSION = "2.0.0";
    constexpr const char* BUILD_DATE = __DATE__;
    constexpr const char* DESCRIPTION = "Molecular modeling data structures for GCMC simulation";
    
    constexpr int TOTAL_COMPONENTS = 8;
    constexpr const char* COMPONENTS[] = {
        "common", "atom", "residue", "molecule", 
        "topology", "montecarlo", "param", "structure"
    };
}

/**
 * @brief Get model module version (following system pattern)
 */
inline std::string getModelVersion() {
    return info::VERSION;
}

/**
 * @brief Get module version and build information
 */
inline std::string get_module_info() {
    std::stringstream ss;
    ss << "Model Module v" << info::VERSION << "\n";
    ss << "Build Date: " << info::BUILD_DATE << "\n";
    ss << "Components: " << info::TOTAL_COMPONENTS << "\n";
    ss << "Description: " << info::DESCRIPTION << "\n";
    return ss.str();
}

/**
 * @brief Module validation and testing utilities
 */
namespace utils {

/**
 * @brief Run comprehensive module validation tests
 */
inline bool validate_module() {
    return factory::run_all_tests();
}

/**
 * @brief Run all factory function tests
 */
inline bool test_all_factories() {
    return factory::test_factory_functions();
}

/**
 * @brief Run basic component tests
 */
inline bool test_basic_components() {
    return factory::run_basic_tests();
}

/**
 * @brief Check if all module components are available
 */
inline bool check_module_integrity() {
    try {
        // Test that all major components can be instantiated
        bool factory_ok = factory::run_basic_tests();
        bool compatibility_ok = factory::test_factory_functions();
        
        return factory_ok && compatibility_ok;
    } catch (const std::exception&) {
        return false;
    }
}

} // namespace utils

/**
 * @brief Backward Compatibility Type Aliases (following system pattern)
 * 
 * These aliases maintain 100% compatibility with existing code that was using
 * individual header files before the modular refactoring.
 */

// === Primary Data Structure Aliases ===
using Atom = atom::Atom;
using Residue = residue::Residue;
using Molecular = molecule::Molecular;
using Structure = structure::Structure;

// === Topology System Aliases ===
using Topology = topology::Topology;
using ForceField = topology::ForceField;

// === Monte Carlo System Aliases ===
using MCState = montecarlo::MCState;

// === Main Parameter Class with Nested Compatibility ===
class Param : public param::Param {
public:
    // Nested type aliases for Python binding compatibility
    using BasicInfo = param::BasicInfo;
    using SpaceInfo = param::SpaceInfo;
    using MCInfo = param::MCParams;
    using EnergyInfo = param::EnergyInfo;
    using FragmentInfo = param::FragmentInfo;
    using BiasInfo = param::BiasInfo;
    using FileInfo = param::FileInfo;
    
    // Inherit all constructors and functionality
    using param::Param::Param;
};

} // namespace model

// Export to parent namespace for convenience (following system pattern)
using model::getModelVersion;

} // namespace pygcmc

/**
 * @brief Module Version Macros
 */
#define PYGCMC_MODEL_VERSION_MAJOR 2
#define PYGCMC_MODEL_VERSION_MINOR 0
#define PYGCMC_MODEL_VERSION_PATCH 0
#define PYGCMC_MODEL_VERSION_STRING "2.0.0"

#endif // PYGCMC_MODEL_UTILS_FACTORY_HPP 