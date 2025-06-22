#pragma once

#ifndef PYGCMC_MODEL_UTILS_VALIDATION_HPP
#define PYGCMC_MODEL_UTILS_VALIDATION_HPP

#include "../atom/AtomMain.hpp"
#include "../residue/ResidueMain.hpp"
#include "../molecule/MolecularMain.hpp"
#include <memory>
#include <string>
#include <sstream>

namespace pygcmc {
namespace model {
namespace validation {

/**
 * @brief Validate a complete molecular system
 * @param system The molecular system to validate
 * @return true if the system is valid, false otherwise
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
 * @brief Validate an atom
 * @param atom The atom to validate
 * @return true if the atom is valid, false otherwise
 */
inline bool validate_atom(const atom::Atom& atom) {
    return atom.is_valid();
}

/**
 * @brief Validate a residue
 * @param residue The residue to validate
 * @return true if the residue is valid, false otherwise
 */
inline bool validate_residue(const residue::Residue& residue) {
    if (!residue.is_valid()) return false;
    
    // Check all atoms in the residue
    for (const auto& atom : residue.get_atoms()) {
        if (!atom || !validate_atom(*atom)) return false;
    }
    
    return true;
}

/**
 * @brief Get system statistics
 * @param system The molecular system to analyze
 * @return A formatted string with system statistics
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
 * @brief Get residue statistics
 * @param residue The residue to analyze
 * @return A formatted string with residue statistics
 */
inline std::string get_residue_summary(const residue::Residue& residue) {
    std::stringstream ss;
    ss << "Residue Summary: " << residue.get_resname() << " (ID: " << residue.get_ires() << ")\n";
    ss << "  Atoms: " << residue.get_atoms().size() << "\n";
    ss << "  Segment: " << residue.get_segid() << "\n";
    
    // Calculate mass and charge
    double total_mass = 0.0;
    double total_charge = 0.0;
    for (const auto& atom : residue.get_atoms()) {
        if (atom) {
            total_mass += atom->get_mass();
            total_charge += atom->get_charge();
        }
    }
    
    ss << "  Total Mass: " << total_mass << " amu\n";
    ss << "  Total Charge: " << total_charge << " e\n";
    return ss.str();
}

/**
 * @brief Check for common molecular system issues
 * @param system The molecular system to check
 * @return A vector of warning/error messages
 */
inline std::vector<std::string> check_system_issues(const molecule::Molecular& system) {
    std::vector<std::string> issues;
    
    // Check for empty system
    if (system.get_atoms().empty()) {
        issues.push_back("Warning: System contains no atoms");
    }
    
    if (system.get_residues().empty()) {
        issues.push_back("Warning: System contains no residues");
    }
    
    // Check for atoms with zero mass
    for (const auto& atom : system.get_atoms()) {
        if (atom && atom->get_mass() <= 0.0) {
            issues.push_back("Error: Atom with zero or negative mass found");
            break;
        }
    }
    
    // Check for extremely large charges
    for (const auto& atom : system.get_atoms()) {
        if (atom && std::abs(atom->get_charge()) > 5.0) {
            issues.push_back("Warning: Atom with unusually large charge found");
            break;
        }
    }
    
    return issues;
}

/**
 * @brief Comprehensive system validation with detailed reporting
 * @param system The molecular system to validate
 * @return A pair containing validation result and detailed report
 */
inline std::pair<bool, std::string> validate_system_detailed(const molecule::Molecular& system) {
    std::stringstream report;
    bool is_valid = true;
    
    report << "=== Molecular System Validation Report ===\n";
    
    // Basic validation
    if (!validate_molecular_system(system)) {
        is_valid = false;
        report << "FAILED: Basic system validation\n";
    } else {
        report << "PASSED: Basic system validation\n";
    }
    
    // Check for issues
    auto issues = check_system_issues(system);
    if (!issues.empty()) {
        report << "\nIssues found:\n";
        for (const auto& issue : issues) {
            report << "  - " << issue << "\n";
            if (issue.find("Error:") != std::string::npos) {
                is_valid = false;
            }
        }
    } else {
        report << "No issues found\n";
    }
    
    // System summary
    report << "\n" << get_system_summary(system);
    
    return {is_valid, report.str()};
}

} // namespace validation
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_UTILS_VALIDATION_HPP 