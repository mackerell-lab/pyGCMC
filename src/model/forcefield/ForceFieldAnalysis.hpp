#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_ANALYSIS_HPP
#define PYGCMC_MODEL_FORCEFIELD_ANALYSIS_HPP

#include "ForceFieldParams.hpp"
#include "ForceFieldInterface.hpp"
#include <sstream>
#include <iomanip>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Force field analysis and statistics class
 * 
 * This class provides analysis and statistics functionality for force fields.
 * It generates summaries, detailed statistics, and performs various analytical
 * operations on force field data.
 */
class ForceFieldAnalyzer {
public:
    /**
     * @brief Constructor
     * @param storage Reference to the force field storage container
     */
    explicit ForceFieldAnalyzer(const ForceFieldStorage& storage) : storage_(storage) {}
    ~ForceFieldAnalyzer() = default;

    // === Statistics Methods ===

    /**
     * @brief Get comprehensive statistics about the force field
     * @return ForceFieldStats structure with counts
     */
    ForceFieldStats get_statistics() const {
        return storage_.get_statistics();
    }

    /**
     * @brief Get a human-readable summary of the force field
     * @return String summary with basic statistics
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
     * @return Detailed string with all parameter listings
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

        // Bonded parameters summary
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

    // === Analysis Methods ===

    /**
     * @brief Find all atom types that have masses but no LJ parameters
     * @return Set of atom types with missing LJ parameters
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
     * @return Set of atom types with missing masses
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
     * @return Set of all atom types found in any parameter set
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

        // From dihedral parameters
        for (const auto& pair : storage_.dihedral_params) {
            types.insert(std::get<0>(pair.first));
            types.insert(std::get<1>(pair.first));
            types.insert(std::get<2>(pair.first));
            types.insert(std::get<3>(pair.first));
        }

        // From improper parameters
        for (const auto& pair : storage_.improper_params) {
            types.insert(std::get<0>(pair.first));
            types.insert(std::get<1>(pair.first));
            types.insert(std::get<2>(pair.first));
            types.insert(std::get<3>(pair.first));
        }

        return types;
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

private:
    const ForceFieldStorage& storage_;  ///< Reference to parameter storage
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_ANALYSIS_HPP 