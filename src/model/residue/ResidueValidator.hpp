#pragma once

#ifndef PYGCMC_MODEL_RESIDUE_VALIDATOR_HPP
#define PYGCMC_MODEL_RESIDUE_VALIDATOR_HPP

#include "ResidueComposite.hpp"
#include "../common/ModelUtils.hpp"
#include <string>
#include <vector>
#include <set>
#include <algorithm>

namespace pygcmc {
namespace model {
namespace residue {

/**
 * @brief Validator class for residue integrity checking
 * Provides comprehensive validation of residue structure and atom composition
 */
class ResidueValidator {
public:
    /**
     * @brief Validation error codes
     */
    enum class ValidationError {
        NONE,
        EMPTY_RESIDUE_NAME,
        INVALID_RESIDUE_NUMBER,
        EMPTY_SEGMENT_ID,
        NO_ATOMS,
        INVALID_ATOM,
        DUPLICATE_ATOM_TYPES,
        INCONSISTENT_RESIDUE_INFO,
        MISSING_BACKBONE_ATOMS,
        INVALID_GEOMETRY,
        CHARGE_IMBALANCE,
        MASS_IMBALANCE
    };

    /**
     * @brief Check if a residue is valid
     */
    static bool is_valid(const ResidueComposite& residue) {
        return get_validation_errors(residue).empty();
    }

    /**
     * @brief Get detailed validation errors
     */
    static std::vector<ValidationError> get_validation_errors(const ResidueComposite& residue) {
        std::vector<ValidationError> errors;

        // Basic field validation
        if (residue.get_resname().empty()) {
            errors.push_back(ValidationError::EMPTY_RESIDUE_NAME);
        }

        if (residue.get_ires() <= 0) {
            errors.push_back(ValidationError::INVALID_RESIDUE_NUMBER);
        }

        if (residue.get_segid().empty()) {
            errors.push_back(ValidationError::EMPTY_SEGMENT_ID);
        }

        // Atom validation
        const auto& atoms = residue.get_atoms();
        if (atoms.empty()) {
            errors.push_back(ValidationError::NO_ATOMS);
            return errors; // No need to check further if no atoms
        }

        // Check each atom
        for (const auto& atom : atoms) {
            if (!atom || !atom->is_valid()) {
                errors.push_back(ValidationError::INVALID_ATOM);
                break; // One invalid atom is enough
            }
        }

        // Check for duplicate atom types
        if (has_duplicate_atom_types(residue)) {
            errors.push_back(ValidationError::DUPLICATE_ATOM_TYPES);
        }

        // Check consistency between residue and atom info
        if (!has_consistent_atom_info(residue)) {
            errors.push_back(ValidationError::INCONSISTENT_RESIDUE_INFO);
        }

        // Check backbone atoms for protein residues
        if (is_protein_residue(residue.get_resname()) && !has_backbone_atoms(residue)) {
            errors.push_back(ValidationError::MISSING_BACKBONE_ATOMS);
        }

        return errors;
    }

    /**
     * @brief Get validation error message
     */
    static std::string get_error_message(ValidationError error) {
        switch (error) {
            case ValidationError::NONE:
                return "No errors";
            case ValidationError::EMPTY_RESIDUE_NAME:
                return "Empty residue name";
            case ValidationError::INVALID_RESIDUE_NUMBER:
                return "Invalid residue number";
            case ValidationError::EMPTY_SEGMENT_ID:
                return "Empty segment ID";
            case ValidationError::NO_ATOMS:
                return "No atoms in residue";
            case ValidationError::INVALID_ATOM:
                return "Contains invalid atom";
            case ValidationError::DUPLICATE_ATOM_TYPES:
                return "Duplicate atom types";
            case ValidationError::INCONSISTENT_RESIDUE_INFO:
                return "Inconsistent residue information";
            case ValidationError::MISSING_BACKBONE_ATOMS:
                return "Missing backbone atoms";
            case ValidationError::INVALID_GEOMETRY:
                return "Invalid geometry";
            case ValidationError::CHARGE_IMBALANCE:
                return "Charge imbalance";
            case ValidationError::MASS_IMBALANCE:
                return "Mass imbalance";
            default:
                return "Unknown error";
        }
    }

    /**
     * @brief Get comprehensive validation report
     */
    static std::string get_validation_report(const ResidueComposite& residue) {
        auto errors = get_validation_errors(residue);
        if (errors.empty()) {
            return "Residue is valid";
        }

        std::string report = "Validation errors:\n";
        for (const auto& error : errors) {
            report += "- " + get_error_message(error) + "\n";
        }
        return report;
    }

    /**
     * @brief Check if residue has duplicate atom types
     */
    static bool has_duplicate_atom_types(const ResidueComposite& residue) {
        const auto& atoms = residue.get_atoms();
        std::set<std::string> atom_types;
        
        for (const auto& atom : atoms) {
            if (!atom) continue;
            const std::string& type = atom->get_type();
            if (atom_types.find(type) != atom_types.end()) {
                return true; // Duplicate found
            }
            atom_types.insert(type);
        }
        return false;
    }

    /**
     * @brief Check if all atoms have consistent residue information
     */
    static bool has_consistent_atom_info(const ResidueComposite& residue) {
        const auto& atoms = residue.get_atoms();
        const std::string& resname = residue.get_resname();
        int ires = residue.get_ires();
        const std::string& segid = residue.get_segid();
        int iseg = residue.get_iseg();

        return std::all_of(atoms.begin(), atoms.end(),
            [&](const std::shared_ptr<atom::Atom>& atom) {
                if (!atom) return false;
                return atom->get_resname() == resname &&
                       atom->get_ires() == ires &&
                       atom->get_segid() == segid &&
                       atom->get_iseg() == iseg;
            });
    }

    /**
     * @brief Check if residue is a protein residue
     */
    static bool is_protein_residue(const std::string& resname) {
        static const std::set<std::string> protein_residues = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", 
            "HIS", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", "MET", 
            "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
        };
        return protein_residues.find(resname) != protein_residues.end();
    }

    /**
     * @brief Check if residue is a nucleic acid residue
     */
    static bool is_nucleic_acid_residue(const std::string& resname) {
        static const std::set<std::string> nucleic_residues = {
            "ADE", "GUA", "CYT", "THY", "URA",  // DNA/RNA bases
            "A", "G", "C", "T", "U"             // Single letter codes
        };
        return nucleic_residues.find(resname) != nucleic_residues.end();
    }

    /**
     * @brief Check if protein residue has backbone atoms
     */
    static bool has_backbone_atoms(const ResidueComposite& residue) {
        static const std::vector<std::string> backbone_atoms = {"N", "CA", "C", "O"};
        
        for (const std::string& atom_type : backbone_atoms) {
            if (!residue.has_atom_type(atom_type)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Check if nucleic acid residue has backbone atoms
     */
    static bool has_nucleic_backbone_atoms(const ResidueComposite& residue) {
        static const std::vector<std::string> nucleic_backbone = {"P", "O5'", "C5'", "C4'", "C3'", "O3'"};
        
        for (const std::string& atom_type : nucleic_backbone) {
            if (!residue.has_atom_type(atom_type)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Calculate total charge of residue
     */
    static double calculate_total_charge(const ResidueComposite& residue) {
        double total_charge = 0.0;
        const auto& atoms = residue.get_atoms();
        
        for (const auto& atom : atoms) {
            if (atom) {
                total_charge += atom->get_charge();
            }
        }
        return total_charge;
    }

    /**
     * @brief Calculate total mass of residue
     */
    static double calculate_total_mass(const ResidueComposite& residue) {
        double total_mass = 0.0;
        const auto& atoms = residue.get_atoms();
        
        for (const auto& atom : atoms) {
            if (atom) {
                total_mass += atom->get_mass();
            }
        }
        return total_mass;
    }

    /**
     * @brief Check if residue charge is within expected range
     */
    static bool is_charge_reasonable(const ResidueComposite& residue, double tolerance = 0.01) {
        double charge = calculate_total_charge(residue);
        
        // For most residues, charge should be close to integer
        double rounded_charge = std::round(charge);
        return std::abs(charge - rounded_charge) < tolerance;
    }

    /**
     * @brief Check if residue mass is within expected range
     */
    static bool is_mass_reasonable(const ResidueComposite& residue) {
        double mass = calculate_total_mass(residue);
        
        // Basic sanity check: mass should be positive and reasonable
        return mass > 10.0 && mass < 1000.0; // Reasonable range for typical residues
    }

    /**
     * @brief Check bond lengths for geometric consistency
     */
    static bool has_reasonable_geometry(const ResidueComposite& residue, double max_bond_length = 5.0) {
        const auto& atoms = residue.get_atoms();
        
        // Check distances between all pairs of atoms
        for (size_t i = 0; i < atoms.size(); ++i) {
            for (size_t j = i + 1; j < atoms.size(); ++j) {
                if (!atoms[i] || !atoms[j]) continue;
                
                double distance = atoms[i]->distance_to(*atoms[j]);
                if (distance > max_bond_length) {
                    // Atoms are too far apart, might indicate problems
                    continue;
                }
                if (distance < 0.5) {
                    // Atoms are too close, indicates overlap
                    return false;
                }
            }
        }
        return true;
    }
};

} // namespace residue
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_RESIDUE_VALIDATOR_HPP 