#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_UTILS_HPP
#define PYGCMC_MODEL_FORCEFIELD_UTILS_HPP

#include "ForceFieldParams.hpp"
#include <sstream>
#include <iomanip>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Force field utility class providing analysis and validation tools
 * 
 * This class provides utility functions for force field analysis, validation,
 * and statistics generation. It operates on a ForceFieldStorage instance
 * to perform various analytical operations.
 */
class ForceFieldUtils {
public:
    /**
     * @brief Constructor
     * @param storage Reference to the force field storage container
     */
    explicit ForceFieldUtils(const ForceFieldStorage& storage) : storage_(storage) {}
    ~ForceFieldUtils() = default;

    // === Validation Methods ===

    /**
     * @brief Validate the force field completeness
     * @return true if force field has basic required parameters
     */
    bool validate_force_field() const {
        // Basic validation: check if we have any atom masses and LJ parameters
        return !storage_.atom_masses.empty() && !storage_.lj_params.empty();
    }

    /**
     * @brief Check completeness for a given set of atom types
     * @param atom_types Set of atom types to check
     * @return CompletenessResult with detailed analysis
     */
    CompletenessResult check_completeness(const std::set<std::string>& atom_types) const {
        CompletenessResult result;
        result.is_complete = true;

        // Check atom masses
        for (const auto& type : atom_types) {
            if (storage_.atom_masses.find(type) == storage_.atom_masses.end()) {
                result.missing_atom_masses.insert(type);
                result.is_complete = false;
            }
        }

        // Check LJ parameters
        for (const auto& type : atom_types) {
            if (storage_.lj_params.find(type) == storage_.lj_params.end()) {
                result.missing_lj_params.insert(type);
                result.is_complete = false;
            }
        }

        // Generate summary
        std::ostringstream oss;
        if (result.is_complete) {
            oss << "Force field is complete for all " << atom_types.size() << " atom types.";
        } else {
            oss << "Force field is incomplete:\n";
            if (!result.missing_atom_masses.empty()) {
                oss << "  Missing atom masses: " << result.missing_atom_masses.size() << " types\n";
            }
            if (!result.missing_lj_params.empty()) {
                oss << "  Missing LJ parameters: " << result.missing_lj_params.size() << " types\n";
            }
        }
        result.summary = oss.str();

        return result;
    }

    // === Statistics Methods ===

    /**
     * @brief Get comprehensive statistics about the force field
     */
    ForceFieldStats get_statistics() const {
        return storage_.get_statistics();
    }

    /**
     * @brief Get a human-readable summary of the force field
     */
    std::string get_summary() const {
        auto stats = get_statistics();
        std::ostringstream oss;
        oss << "Force Field Summary:\n";
        oss << "  Atom types: " << stats.num_atom_types << "\n";
        oss << "  LJ parameters: " << stats.num_lj_params << "\n";
        oss << "  NBFIX entries: " << stats.num_nbfix << "\n";
        oss << "  Bond types: " << stats.num_bond_types << "\n";
        oss << "  Angle types: " << stats.num_angle_types << "\n";
        oss << "  Dihedral types: " << stats.num_dihedral_types << "\n";
        oss << "  Improper types: " << stats.num_improper_types;
        return oss.str();
    }

    /**
     * @brief Get detailed statistics with parameter information
     */
    std::string get_detailed_statistics() const {
        std::ostringstream oss;
        auto stats = get_statistics();

        oss << "=== Detailed Force Field Statistics ===\n\n";

        // Atom masses
        oss << "Atom Masses (" << stats.num_atom_types << " types):\n";
        for (const auto& pair : storage_.atom_masses) {
            oss << "  " << std::setw(8) << pair.first << ": " 
                << std::setw(8) << std::fixed << std::setprecision(3) << pair.second << " amu\n";
        }
        oss << "\n";

        // LJ parameters
        oss << "Lennard-Jones Parameters (" << stats.num_lj_params << " types):\n";
        for (const auto& pair : storage_.lj_params) {
            oss << "  " << std::setw(8) << pair.first << ": " 
                << "eps=" << std::setw(8) << std::fixed << std::setprecision(4) << pair.second.epsilon
                << " rmin/2=" << std::setw(8) << std::fixed << std::setprecision(4) << pair.second.rmin_half << "\n";
        }
        oss << "\n";

        // NBFIX entries
        if (stats.num_nbfix > 0) {
            oss << "NBFIX Entries (" << stats.num_nbfix << "):\n";
            for (const auto& pair : storage_.nbfix) {
                oss << "  " << std::setw(8) << pair.first.first << "-" << std::setw(8) << pair.first.second << ": "
                    << "eps=" << std::setw(8) << std::fixed << std::setprecision(4) << pair.second.epsilon
                    << " rmin=" << std::setw(8) << std::fixed << std::setprecision(4) << pair.second.rmin << "\n";
            }
            oss << "\n";
        }

        // Nonbonded parameters
        oss << "Nonbonded Parameters:\n";
        const auto& nb = storage_.nonbonded_params;
        oss << "  nbxmod: " << nb.nbxmod << "\n";
        oss << "  cutnb: " << nb.cutnb << " Å\n";
        oss << "  ctofnb: " << nb.ctofnb << " Å\n";
        oss << "  ctonnb: " << nb.ctonnb << " Å\n";
        oss << "  eps: " << nb.eps << "\n";
        oss << "  e14fac: " << nb.e14fac << "\n";
        oss << "  Flags: ";
        if (nb.cdiel) oss << "cdiel ";
        if (nb.fshift) oss << "fshift ";
        if (nb.vatom) oss << "vatom ";
        if (nb.vdistance) oss << "vdistance ";
        if (nb.vfswitch) oss << "vfswitch ";
        oss << "\n\n";

        // Bond parameters summary
        if (stats.num_bond_types > 0) {
            oss << "Bond Parameters: " << stats.num_bond_types << " types\n";
        }
        if (stats.num_angle_types > 0) {
            oss << "Angle Parameters: " << stats.num_angle_types << " types\n";
        }
        if (stats.num_dihedral_types > 0) {
            oss << "Dihedral Parameters: " << stats.num_dihedral_types << " types\n";
        }
        if (stats.num_improper_types > 0) {
            oss << "Improper Parameters: " << stats.num_improper_types << " types\n";
        }

        return oss.str();
    }

    // === Direct Access Methods (for Python bindings) ===

    const std::map<std::string, double>& get_atom_masses() const { 
        return storage_.atom_masses; 
    }
    
    const std::map<std::string, LJParams>& get_lj_params() const { 
        return storage_.lj_params; 
    }
    
    const std::map<std::pair<std::string, std::string>, NBFIXParams>& get_nbfix() const { 
        return storage_.nbfix; 
    }
    
    const std::map<std::pair<std::string, std::string>, BondParams>& get_bond_params() const { 
        return storage_.bond_params; 
    }
    
    const std::map<std::tuple<std::string, std::string, std::string>, AngleParams>& get_angle_params() const { 
        return storage_.angle_params; 
    }
    
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>>& get_dihedral_params() const { 
        return storage_.dihedral_params; 
    }
    
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams>& get_improper_params() const { 
        return storage_.improper_params; 
    }

    // === Analysis Methods ===

    /**
     * @brief Find all atom types that have masses but no LJ parameters
     */
    std::set<std::string> find_missing_lj_params() const {
        std::set<std::string> missing;
        for (const auto& pair : storage_.atom_masses) {
            if (storage_.lj_params.find(pair.first) == storage_.lj_params.end()) {
                missing.insert(pair.first);
            }
        }
        return missing;
    }

    /**
     * @brief Find all atom types that have LJ parameters but no masses
     */
    std::set<std::string> find_missing_masses() const {
        std::set<std::string> missing;
        for (const auto& pair : storage_.lj_params) {
            if (storage_.atom_masses.find(pair.first) == storage_.atom_masses.end()) {
                missing.insert(pair.first);
            }
        }
        return missing;
    }

    /**
     * @brief Get all unique atom types referenced in the force field
     */
    std::set<std::string> get_all_referenced_types() const {
        std::set<std::string> types;

        // From atom masses
        for (const auto& pair : storage_.atom_masses) {
            types.insert(pair.first);
        }

        // From LJ parameters
        for (const auto& pair : storage_.lj_params) {
            types.insert(pair.first);
        }

        // From NBFIX
        for (const auto& pair : storage_.nbfix) {
            types.insert(pair.first.first);
            types.insert(pair.first.second);
        }

        // From bond parameters
        for (const auto& pair : storage_.bond_params) {
            types.insert(pair.first.first);
            types.insert(pair.first.second);
        }

        // From angle parameters
        for (const auto& pair : storage_.angle_params) {
            types.insert(std::get<0>(pair.first));
            types.insert(std::get<1>(pair.first));
            types.insert(std::get<2>(pair.first));
        }

        return types;
    }

    /**
     * @brief Check if force field parameters are consistent
     */
    std::vector<std::string> validate_consistency() const {
        std::vector<std::string> issues;

        // Check for missing masses
        auto missing_masses = find_missing_masses();
        if (!missing_masses.empty()) {
            std::ostringstream oss;
            oss << "Missing atom masses for " << missing_masses.size() << " types with LJ parameters";
            issues.push_back(oss.str());
        }

        // Check for missing LJ parameters
        auto missing_lj = find_missing_lj_params();
        if (!missing_lj.empty()) {
            std::ostringstream oss;
            oss << "Missing LJ parameters for " << missing_lj.size() << " types with masses";
            issues.push_back(oss.str());
        }

        // Check nonbonded parameter consistency
        const auto& nb = storage_.nonbonded_params;
        if (nb.cutnb <= nb.ctofnb) {
            issues.push_back("cutnb should be greater than ctofnb");
        }
        if (nb.ctofnb <= nb.ctonnb) {
            issues.push_back("ctofnb should be greater than ctonnb");
        }

        return issues;
    }

private:
    const ForceFieldStorage& storage_;
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_UTILS_HPP 