#pragma once

#include "ResidueCore.hpp"
#include "../atom/AtomMain.hpp"
#include <vector>
#include <memory>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <stdexcept>
#include <cmath>

namespace pygcmc {
namespace model {
namespace residue {

/**
 * @brief Residue class following CHARMM naming conventions
 */
class Residue {
public:
    // Constructors
    Residue() = default;

    Residue(const std::string& resname, int ires, 
            const std::string& segid = "", int iseg = 0,
            char chain = ' ', char inscode = ' ') {
        data_.resname = resname;
        data_.ires = ires;
        data_.segid = segid;
        data_.iseg = iseg;
        data_.chain = chain;
        data_.inscode = inscode;
    }

    // CHARMM standard getters
    const std::string& get_resname() const noexcept { return data_.resname; }
    int get_ires() const noexcept { return data_.ires; }
    const std::string& get_segid() const noexcept { return data_.segid; }
    int get_iseg() const noexcept { return data_.iseg; }
    char get_chain() const noexcept { return data_.chain; }
    char get_inscode() const noexcept { return data_.inscode; }
    
    // CHARMM standard setters
    void set_resname(const std::string& name) noexcept { data_.resname = name; }
    
    // Atom management
    void add_atom(const atom::Atom& atom) {
        // Verify atom belongs to this residue
        if (atom.get_resname() != data_.resname || atom.get_ires() != data_.ires ||
            atom.get_segid() != data_.segid || atom.get_iseg() != data_.iseg) {
            throw std::invalid_argument("Atom does not belong to this residue");
        }
        atoms_.push_back(std::make_shared<atom::Atom>(atom));
    }

    void add_atom(std::shared_ptr<atom::Atom> atom) {
        if (!atom) return;
        if (atom->get_resname() != data_.resname || atom->get_ires() != data_.ires ||
            atom->get_segid() != data_.segid || atom->get_iseg() != data_.iseg) {
            throw std::invalid_argument("Atom does not belong to this residue");
        }
        atoms_.push_back(atom);
    }

    const std::vector<std::shared_ptr<atom::Atom>>& get_atoms() const noexcept { 
        return atoms_; 
    }

    std::shared_ptr<atom::Atom> find_atom(const std::string& type) const {
        auto it = std::find_if(atoms_.begin(), atoms_.end(),
            [&type](const std::shared_ptr<atom::Atom>& atom) {
                return atom && atom->get_type() == type;
            });
        return (it != atoms_.end()) ? *it : nullptr;
    }

    // Utility methods
    size_t atom_count() const noexcept {
        return atoms_.size();
    }

    bool is_valid() const {
        return !data_.resname.empty() && data_.ires > 0 && !data_.segid.empty() &&
               std::all_of(atoms_.begin(), atoms_.end(),
                          [](const std::shared_ptr<atom::Atom>& atom) {
                              return atom && atom->is_valid();
                          });
    }

    // Center of mass calculation and storage
    void calculate_center_of_mass() {
        data_.com = {0.0, 0.0, 0.0};
        double totalMass = 0.0;

        for (const auto& atom : atoms_) {
            if (!atom) continue;
            double mass = atom->get_mass();
            const auto& coor = atom->get_coor();
            for (int i = 0; i < 3; ++i) {
                data_.com[i] += mass * coor[i];
            }
            totalMass += mass;
        }

        if (totalMass > 0.0) {
            for (double& x : data_.com) x /= totalMass;
        }
    }

    const std::array<double, 3>& get_center_of_mass() const noexcept {
        return data_.com;
    }

    // Selection methods for CHARMM compatibility
    bool has_atom_type(const std::string& type) const {
        return std::any_of(atoms_.begin(), atoms_.end(),
            [&type](const std::shared_ptr<atom::Atom>& atom) {
                return atom && atom->get_type() == type;
            });
    }

    std::vector<std::shared_ptr<atom::Atom>> select_atoms(
        const std::function<bool(const atom::Atom&)>& predicate) const {
        std::vector<std::shared_ptr<atom::Atom>> selected;
        for (const auto& atom : atoms_) {
            if (atom && predicate(*atom)) {
                selected.push_back(atom);
            }
        }
        return selected;
    }

    // Secondary structure getters/setters
    SecondaryStructure get_secondary_structure() const noexcept { return data_.secStruct; }
    const SheetStrand& get_sheet_info() const noexcept { return data_.sheetInfo; }
    const SSBond& get_ssbond() const noexcept { return data_.ssbond; }

    void set_secondary_structure(SecondaryStructure ss) { data_.secStruct = ss; }
    void set_sheet_info(const SheetStrand& si) { data_.sheetInfo = si; }
    void set_ssbond(const SSBond& sb) { data_.ssbond = sb; }

    // PDB format utilities
    std::string get_residue_id() const {
        // Combine residue number and insertion code (e.g., "153A")
        if (data_.inscode == ' ') {
            return std::to_string(data_.ires);
        }
        return std::to_string(data_.ires) + data_.inscode;
    }

    void set_residue_id(const std::string& resid) {
        // Parse residue ID (e.g., "153A" -> ires=153, inscode='A')
        size_t numLen = 0;
        try {
            data_.ires = std::stoi(resid, &numLen);
        } catch (const std::exception&) {
            throw std::invalid_argument("Invalid residue ID format");
        }
        
        if (numLen < resid.length()) {
            data_.inscode = resid[numLen];
        } else {
            data_.inscode = ' ';
        }
    }

    // Enhanced atom lookup methods
    std::shared_ptr<atom::Atom> find_atom_by_pdb_name(const std::string& pdbName) const {
        // Find atom by PDB formatted name
        auto it = std::find_if(atoms_.begin(), atoms_.end(),
            [&pdbName](const std::shared_ptr<atom::Atom>& atom) {
                return atom && atom->get_formatted_atom_name() == pdbName;
            });
        return (it != atoms_.end()) ? *it : nullptr;
    }

    // CHARMM-style atom range
    std::pair<size_t, size_t> get_atom_range() const {
        return {0, atoms_.size()};  // Equivalent to IBASE(IRES) to IBASE(IRES+1)
    }

    // HETATM support
    bool is_hetatm() const noexcept { return data_.hetatm; }
    void set_hetatm(bool het) noexcept { data_.hetatm = het; }





private:
    void update_atom_map() {
        atomMap_.clear();
        for (const auto& atom : atoms_) {
            if (atom) {
                // Store both raw and PDB-formatted names
                atomMap_[atom->get_type()] = atom;
                atomMap_[atom->get_formatted_atom_name()] = atom;
            }
        }
    }

    // Core data
    ResidueData data_;

    // Atom storage
    std::vector<std::shared_ptr<atom::Atom>> atoms_;  ///< Atoms in residue
    std::unordered_map<std::string, std::shared_ptr<atom::Atom>> atomMap_;  ///< Quick atom lookup by type
};

/**
 * @brief Validate a residue (from original codebase)
 * @param residue The residue to validate
 * @return true if the residue is valid, false otherwise
 */
inline bool validate_residue(const Residue& residue) {
    if (!residue.is_valid()) return false;
    
    // Check all atoms in the residue
    for (const auto& atom : residue.get_atoms()) {
        if (!atom || !atom->is_valid()) return false;
    }
    
    return true;
}

} // namespace residue
} // namespace model
} // namespace pygcmc