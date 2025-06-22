#pragma once

#ifndef PYGCMC_MODEL_UTILS_FACTORY_HPP
#define PYGCMC_MODEL_UTILS_FACTORY_HPP

#include "../atom/AtomMain.hpp"
#include "../residue/ResidueMain.hpp"
#include "../molecule/MolecularMain.hpp"
#include <memory>
#include <string>
#include <sstream>

namespace utils {

/**
 * @brief Create a water molecule
 */
inline std::shared_ptr<pygcmc::model::residue::Residue> create_water_molecule(int residue_id, const std::string& segment = "SOLV") {
    auto water = std::make_shared<pygcmc::model::residue::Residue>("SOL", residue_id, segment);
    
    // Add atoms: O, H1, H2 with correct constructor parameters
    // AtomCore(int bynu, type, resname, ires, segid, iseg, x, y, z, wmain, mass, charge, chem, hetatm)
    auto oxygen = std::make_shared<pygcmc::model::atom::Atom>(1, "OT", "SOL", residue_id, segment, 0,
                                               0.0, 0.0, 0.0, 1.0, 15.999, -0.834, "O", false);
    auto h1 = std::make_shared<pygcmc::model::atom::Atom>(2, "HT", "SOL", residue_id, segment, 0,
                                           0.757, 0.586, 0.0, 1.0, 1.008, 0.417, "H", false);
    auto h2 = std::make_shared<pygcmc::model::atom::Atom>(3, "HT", "SOL", residue_id, segment, 0,
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
inline std::shared_ptr<pygcmc::model::residue::Residue> create_alanine_residue(int residue_id, const std::string& segment = "PROT") {
    auto ala = std::make_shared<pygcmc::model::residue::Residue>("ALA", residue_id, segment);
    
    // Add backbone atoms with correct constructor parameters
    auto n = std::make_shared<pygcmc::model::atom::Atom>(1, "NH1", "ALA", residue_id, segment, 0,
                                          0.0, 0.0, 0.0, 1.0, 14.007, -0.47, "N", false);
    auto ca = std::make_shared<pygcmc::model::atom::Atom>(2, "CT1", "ALA", residue_id, segment, 0,
                                           1.458, 0.0, 0.0, 1.0, 12.01, 0.07, "C", false);
    auto c = std::make_shared<pygcmc::model::atom::Atom>(3, "C", "ALA", residue_id, segment, 0,
                                          2.009, 1.421, 0.0, 1.0, 12.01, 0.51, "C", false);
    auto o = std::make_shared<pygcmc::model::atom::Atom>(4, "O", "ALA", residue_id, segment, 0,
                                          1.239, 2.364, 0.0, 1.0, 15.999, -0.51, "O", false);
    auto cb = std::make_shared<pygcmc::model::atom::Atom>(5, "CT3", "ALA", residue_id, segment, 0,
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
inline bool validate_molecular_system(const pygcmc::model::molecule::Molecular& system) {
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
inline std::string get_system_summary(const pygcmc::model::molecule::Molecular& system) {
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
        auto atom = std::make_shared<pygcmc::model::atom::Atom>(1, "CT1", "ALA", 1, "PROT", 0, 
                                                  0.0, 0.0, 0.0, 1.0, 12.01, 0.07, "C", false);
        if (!atom->is_valid()) return false;
        
        // Test residue creation
        auto residue = create_alanine_residue(1);
        if (!residue->is_valid()) return false;
        
        // Test molecular system
        auto molecular = std::make_shared<pygcmc::model::molecule::Molecular>();
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
        auto molecular = std::make_shared<pygcmc::model::molecule::Molecular>();
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

} // namespace utils

#endif // PYGCMC_MODEL_UTILS_FACTORY_HPP 