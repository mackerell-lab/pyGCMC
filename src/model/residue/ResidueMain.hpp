#pragma once

#ifndef PYGCMC_MODEL_RESIDUE_MAIN_HPP
#define PYGCMC_MODEL_RESIDUE_MAIN_HPP

#include "ResidueComposite.hpp"
#include "ResidueValidator.hpp"
#include <sstream>
#include <iomanip>

namespace pygcmc {
namespace model {

/**
 * @brief Extended Residue class with validation and PDB compatibility
 * @details Provides complete residue functionality while maintaining backward compatibility
 */
class Residue : public ResidueComposite, public ICloneable<Residue>, public ISerializable {
public:
    // Sheet strand information structure
    struct SheetStrand {
        int strandNumber = 0;           // Strand number in current sheet
        std::string sheetId;            // Sheet identifier
        int numStrands = 0;             // Number of strands in current sheet
        int sense = 0;                  // 0: first, 1: parallel, -1: antiparallel
        // Hydrogen bond information
        std::string atom1;              // First atom name
        std::string atom2;              // Second atom name
        int resnum1 = 0;                // First residue number
        int resnum2 = 0;                // Second residue number
    };

    // Disulfide bond information structure
    struct SSBond {
        int serialNumber = 0;           // Bond serial number
        std::string partner;            // Partner residue identifier
        char partnerChain = ' ';        // Partner chain identifier
        int partnerResnum = 0;          // Partner residue number
        char partnerInscode = ' ';      // Partner insertion code
        double bondLength = 0.0;        // Length of the disulfide bond
    };

    // Inherit constructors
    using ResidueComposite::ResidueComposite;

    // Default constructor
    Residue() = default;

    // Copy constructor and assignment
    Residue(const Residue&) = default;
    Residue& operator=(const Residue&) = default;
    Residue(Residue&&) = default;
    Residue& operator=(Residue&&) = default;

    // Constructor from ResidueComposite
    explicit Residue(const ResidueComposite& composite) : ResidueComposite(composite) {}

    // Extended getters
    const SheetStrand& get_sheet_info() const noexcept { return sheet_info_; }
    const SSBond& get_ssbond() const noexcept { return ssbond_; }

    // Extended setters
    void set_sheet_info(const SheetStrand& si) { sheet_info_ = si; }
    void set_ssbond(const SSBond& sb) { ssbond_ = sb; }

    // PDB format utilities
    std::string get_residue_id() const {
        if (get_inscode() == ' ') {
            return std::to_string(get_ires());
        }
        return std::to_string(get_ires()) + get_inscode();
    }

    void set_residue_id(const std::string& resid) {
        size_t numLen = 0;
        try {
            int ires = std::stoi(resid, &numLen);
            set_ires(ires);
        } catch (const std::exception&) {
            throw std::invalid_argument("Invalid residue ID format");
        }
        
        if (numLen < resid.length()) {
            set_inscode(resid[numLen]);
        } else {
            set_inscode(' ');
        }
    }

    // Enhanced atom lookup methods
    std::shared_ptr<Atom> find_atom_by_pdb_name(const std::string& pdb_name) const {
        auto atoms = get_atoms();
        auto it = std::find_if(atoms.begin(), atoms.end(),
            [&pdb_name](const std::shared_ptr<Atom>& atom) {
                return atom && atom->get_formatted_atom_name() == pdb_name;
            });
        return (it != atoms.end()) ? *it : nullptr;
    }

    // CHARMM-style atom range
    std::pair<size_t, size_t> get_atom_range() const {
        return {0, atom_count()};  // Equivalent to IBASE(IRES) to IBASE(IRES+1)
    }

    // Validation with detailed reporting
    ResidueValidator::ValidationResult validate_with_details() const {
        ResidueValidator validator;
        return validator.validate_residue(*this);
    }

    bool is_complete() const {
        ResidueValidator validator;
        auto result = validator.check_completeness(*this);
        return result.is_valid;
    }

    std::vector<std::string> get_missing_atoms() const {
        ResidueValidator validator;
        return validator.get_missing_atoms(*this);
    }

    std::vector<std::string> get_unexpected_atoms() const {
        ResidueValidator validator;
        return validator.get_unexpected_atoms(*this);
    }

    // Distance calculations
    double distance_to(const Residue& other) const {
        const auto& com1 = get_center_of_mass();
        const auto& com2 = other.get_center_of_mass();
        return utils::math::distance(com1, com2);
    }

    double min_atom_distance_to(const Residue& other) const {
        double min_dist = std::numeric_limits<double>::max();
        
        for (const auto& atom1 : get_atoms()) {
            if (!atom1) continue;
            for (const auto& atom2 : other.get_atoms()) {
                if (!atom2) continue;
                double dist = atom1->distance_to(*atom2);
                min_dist = std::min(min_dist, dist);
            }
        }
        
        return min_dist;
    }

    // Comparison operators
    bool operator==(const Residue& other) const {
        return get_resname() == other.get_resname() &&
               get_ires() == other.get_ires() &&
               get_segid() == other.get_segid() &&
               get_chain() == other.get_chain() &&
               atom_count() == other.atom_count();
    }

    bool operator!=(const Residue& other) const {
        return !(*this == other);
    }

    // Hash support
    std::size_t hash() const {
        return utils::hash::combine_hash(
            utils::hash::string_hash(get_resname()),
            get_ires(),
            utils::hash::string_hash(get_segid()),
            static_cast<std::size_t>(get_chain()),
            atom_count()
        );
    }

    // ICloneable interface
    std::unique_ptr<Residue> clone() const override {
        return std::make_unique<Residue>(*this);
    }

    // ISerializable interface
    std::string serialize() const override {
        std::ostringstream oss;
        oss << "RESIDUE:" << get_resname() << ":" << get_ires() 
            << ":" << get_segid() << ":" << get_chain() << ":" << atom_count();
        
        // Serialize atoms
        for (const auto& atom : get_atoms()) {
            if (atom) {
                oss << "|" << atom->serialize();
            }
        }
        
        return oss.str();
    }

    bool deserialize(const std::string& data) override {
        std::istringstream iss(data);
        std::string header;
        
        if (!std::getline(iss, header, '|')) return false;
        
        // Parse header: RESIDUE:resname:ires:segid:chain:atom_count
        std::istringstream header_stream(header);
        std::string token;
        
        if (!std::getline(header_stream, token, ':') || token != "RESIDUE") return false;
        
        std::string resname;
        if (!std::getline(header_stream, resname, ':')) return false;
        set_resname(resname);
        
        if (!std::getline(header_stream, token, ':')) return false;
        set_ires(std::stoi(token));
        
        std::string segid;
        if (!std::getline(header_stream, segid, ':')) return false;
        set_segid(segid);
        
        if (!std::getline(header_stream, token, ':')) return false;
        set_chain(token.empty() ? ' ' : token[0]);
        
        if (!std::getline(header_stream, token)) return false;
        size_t expected_atoms = std::stoull(token);
        
        // Deserialize atoms
        clear_atoms();
        std::string atom_data;
        for (size_t i = 0; i < expected_atoms && std::getline(iss, atom_data, '|'); ++i) {
            auto atom = std::make_shared<Atom>();
            if (atom->deserialize(atom_data)) {
                add_atom(atom);
            }
        }
        
        return atom_count() == expected_atoms;
    }

    std::string get_type_name() const override {
        return "Residue";
    }

    // PDB format output
    std::string to_pdb_string() const {
        std::ostringstream oss;
        
        for (const auto& atom : get_atoms()) {
            if (atom) {
                oss << atom->to_pdb_string() << "\n";
            }
        }
        
        return oss.str();
    }

    // Enhanced validation with detailed error messages
    std::string get_validation_error() const override {
        auto result = validate_with_details();
        if (result.is_valid) {
            return "";
        }
        
        std::ostringstream oss;
        for (const auto& error : result.errors) {
            if (!oss.str().empty()) oss << "; ";
            oss << error;
        }
        
        return oss.str();
    }

private:
    SheetStrand sheet_info_;    ///< Sheet strand information
    SSBond ssbond_;             ///< Disulfide bond information
};

} // namespace model
} // namespace pygcmc

// Hash specialization for std::unordered_map support
namespace std {
    template<>
    struct hash<pygcmc::model::Residue> {
        std::size_t operator()(const pygcmc::model::Residue& residue) const {
            return residue.hash();
        }
    };
}

#endif // PYGCMC_MODEL_RESIDUE_MAIN_HPP 