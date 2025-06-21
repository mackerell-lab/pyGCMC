#pragma once

#ifndef PYGCMC_MODEL_RESIDUE_VALIDATOR_HPP
#define PYGCMC_MODEL_RESIDUE_VALIDATOR_HPP

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include "../atom/AtomMain.hpp"
#include "../common/ModelUtils.hpp"

namespace pygcmc {
namespace model {

// Forward declaration
class ResidueComposite;

/**
 * @brief Residue validation and completeness checking utilities
 */
class ResidueValidator {
public:
    // Validation result structure
    struct ValidationResult {
        bool is_valid = true;
        std::vector<std::string> errors;
        std::vector<std::string> warnings;
        
        void add_error(const std::string& error) {
            errors.push_back(error);
            is_valid = false;
        }
        
        void add_warning(const std::string& warning) {
            warnings.push_back(warning);
        }
        
        std::string get_summary() const {
            std::string result = "Validation " + (is_valid ? "PASSED" : "FAILED");
            if (!errors.empty()) {
                result += "\nErrors:";
                for (const auto& error : errors) {
                    result += "\n  - " + error;
                }
            }
            if (!warnings.empty()) {
                result += "\nWarnings:";
                for (const auto& warning : warnings) {
                    result += "\n  - " + warning;
                }
            }
            return result;
        }
    };

    // Standard amino acid definitions
    struct AminoAcidTemplate {
        std::string name;
        std::unordered_set<std::string> required_atoms;
        std::unordered_set<std::string> optional_atoms;
        std::unordered_map<std::string, std::string> atom_elements;
        double expected_mass;
    };

    // Standard nucleotide definitions
    struct NucleotideTemplate {
        std::string name;
        std::unordered_set<std::string> required_atoms;
        std::unordered_set<std::string> optional_atoms;
        double expected_mass;
    };

    ResidueValidator() {
        initialize_templates();
    }

    /**
     * @brief Validate a complete residue
     */
    ValidationResult validate_residue(const ResidueComposite& residue) const;

    /**
     * @brief Check if residue has all required atoms
     */
    ValidationResult check_completeness(const ResidueComposite& residue) const;

    /**
     * @brief Validate atom consistency within residue
     */
    ValidationResult check_atom_consistency(const ResidueComposite& residue) const;

    /**
     * @brief Check for duplicate atoms
     */
    ValidationResult check_duplicates(const ResidueComposite& residue) const;

    /**
     * @brief Validate chemical consistency
     */
    ValidationResult check_chemistry(const ResidueComposite& residue) const;

    /**
     * @brief Check coordinate validity
     */
    ValidationResult check_coordinates(const ResidueComposite& residue) const;

    /**
     * @brief Check mass and charge consistency
     */
    ValidationResult check_mass_charge(const ResidueComposite& residue) const;

    /**
     * @brief Validate residue naming conventions
     */
    ValidationResult check_naming_conventions(const ResidueComposite& residue) const;

    /**
     * @brief Check backbone connectivity (for proteins)
     */
    ValidationResult check_backbone_connectivity(const ResidueComposite& residue) const;

    /**
     * @brief Add custom residue template
     */
    void add_amino_acid_template(const AminoAcidTemplate& template_def) {
        amino_acid_templates_[template_def.name] = template_def;
    }

    void add_nucleotide_template(const NucleotideTemplate& template_def) {
        nucleotide_templates_[template_def.name] = template_def;
    }

    /**
     * @brief Check if residue type is known
     */
    bool is_known_amino_acid(const std::string& resname) const {
        return amino_acid_templates_.find(resname) != amino_acid_templates_.end();
    }

    bool is_known_nucleotide(const std::string& resname) const {
        return nucleotide_templates_.find(resname) != nucleotide_templates_.end();
    }

    /**
     * @brief Get missing atoms for a residue
     */
    std::vector<std::string> get_missing_atoms(const ResidueComposite& residue) const;

    /**
     * @brief Get unexpected atoms for a residue
     */
    std::vector<std::string> get_unexpected_atoms(const ResidueComposite& residue) const;

private:
    std::unordered_map<std::string, AminoAcidTemplate> amino_acid_templates_;
    std::unordered_map<std::string, NucleotideTemplate> nucleotide_templates_;

    void initialize_templates();
    void initialize_amino_acid_templates();
    void initialize_nucleotide_templates();

    ValidationResult validate_against_amino_acid_template(
        const ResidueComposite& residue, 
        const AminoAcidTemplate& template_def) const;

    ValidationResult validate_against_nucleotide_template(
        const ResidueComposite& residue,
        const NucleotideTemplate& template_def) const;

    bool check_bond_distance(const Atom& atom1, const Atom& atom2, 
                           double min_dist, double max_dist) const;
};

// Implementation of validation methods
inline void ResidueValidator::initialize_templates() {
    initialize_amino_acid_templates();
    initialize_nucleotide_templates();
}

inline void ResidueValidator::initialize_amino_acid_templates() {
    // Standard amino acids with required atoms
    amino_acid_templates_["ALA"] = {
        "ALA", {"N", "CA", "C", "O", "CB"}, {"H", "HA", "HB1", "HB2", "HB3"},
        {{"N", "N"}, {"CA", "C"}, {"C", "C"}, {"O", "O"}, {"CB", "C"}}, 89.09
    };
    
    amino_acid_templates_["GLY"] = {
        "GLY", {"N", "CA", "C", "O"}, {"H", "HA2", "HA3"},
        {{"N", "N"}, {"CA", "C"}, {"C", "C"}, {"O", "O"}}, 75.07
    };
    
    amino_acid_templates_["VAL"] = {
        "VAL", {"N", "CA", "C", "O", "CB", "CG1", "CG2"}, 
        {"H", "HA", "HB", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23"},
        {{"N", "N"}, {"CA", "C"}, {"C", "C"}, {"O", "O"}, {"CB", "C"}, {"CG1", "C"}, {"CG2", "C"}}, 117.15
    };
    
    // Add more amino acids as needed...
}

inline void ResidueValidator::initialize_nucleotide_templates() {
    // Standard nucleotides
    nucleotide_templates_["DA"] = {
        "DA", {"P", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'", "N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"},
        {"H5'", "H5''", "H4'", "H3'", "H2'", "H2''", "H1'", "H8", "H61", "H62", "H2"}, 331.22
    };
    
    // Add more nucleotides as needed...
}

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_RESIDUE_VALIDATOR_HPP 