#pragma once

#ifndef PYGCMC_MODEL_RESIDUE_COMPOSITE_HPP
#define PYGCMC_MODEL_RESIDUE_COMPOSITE_HPP

#include <vector>
#include <memory>
#include <string>
#include <array>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include "../atom/AtomMain.hpp"
#include "../common/ModelInterface.hpp"
#include "../common/ModelUtils.hpp"

namespace pygcmc {
namespace model {

/**
 * @brief Core residue functionality for atom composition and management
 */
class ResidueComposite : public IValidatable, public IIdentifiable {
public:
    // Secondary structure types
    enum class SecondaryStructure {
        NONE,
        HELIX_ALPHA_RIGHT,    // Right-handed alpha (default)
        HELIX_OMEGA_RIGHT,    // Right-handed omega
        HELIX_PI_RIGHT,       // Right-handed pi
        HELIX_GAMMA_RIGHT,    // Right-handed gamma
        HELIX_310_RIGHT,      // Right-handed 3/10
        HELIX_ALPHA_LEFT,     // Left-handed alpha
        HELIX_OMEGA_LEFT,     // Left-handed omega
        HELIX_GAMMA_LEFT,     // Left-handed gamma
        HELIX_27_RIBBON,      // 2/7 ribbon/helix
        HELIX_POLYPROLINE,    // Polyproline
        SHEET_PARALLEL,       // Parallel beta sheet
        SHEET_ANTIPARALLEL    // Anti-parallel beta sheet
    };

    // Default constructor
    ResidueComposite() : 
        resname_(""), ires_(0), segid_(""), iseg_(0), igro_(0),
        chain_(' '), inscode_(' '), move_(1), ignore_(0), constrain_(0),
        hetatm_(false), sec_struct_(SecondaryStructure::NONE) {}

    // Parameterized constructor
    ResidueComposite(const std::string& resname, int ires, 
                    const std::string& segid = "", int iseg = 0,
                    char chain = ' ', char inscode = ' ') :
        resname_(resname), ires_(ires), segid_(segid), iseg_(iseg), igro_(0),
        chain_(chain), inscode_(inscode), move_(1), ignore_(0), constrain_(0),
        hetatm_(false), sec_struct_(SecondaryStructure::NONE) {}

    // Basic getters
    const std::string& get_resname() const noexcept { return resname_; }
    int get_ires() const noexcept { return ires_; }
    const std::string& get_segid() const noexcept { return segid_; }
    int get_iseg() const noexcept { return iseg_; }
    char get_chain() const noexcept { return chain_; }
    char get_inscode() const noexcept { return inscode_; }
    SecondaryStructure get_secondary_structure() const noexcept { return sec_struct_; }
    bool is_hetatm() const noexcept { return hetatm_; }

    // Basic setters
    void set_resname(const std::string& name) {
        if (!utils::validate::is_valid_name(name, constants::MAX_RESIDUE_NAME_LENGTH)) {
            throw std::invalid_argument("Invalid residue name");
        }
        resname_ = name;
    }

    void set_ires(int ires) {
        if (ires <= 0) throw std::invalid_argument("Invalid residue number");
        ires_ = ires;
    }

    void set_segid(const std::string& segid) {
        if (segid.length() > constants::MAX_SEGMENT_NAME_LENGTH) {
            throw std::invalid_argument("Segment ID too long");
        }
        segid_ = segid;
    }

    void set_chain(char chain) { chain_ = chain; }
    void set_inscode(char inscode) { inscode_ = inscode; }
    void set_secondary_structure(SecondaryStructure ss) { sec_struct_ = ss; }
    void set_hetatm(bool hetatm) { hetatm_ = hetatm; }

    // Atom management
    void add_atom(const Atom& atom) {
        // Verify atom belongs to this residue
        if (atom.get_resname() != resname_ || atom.get_ires() != ires_ ||
            atom.get_segid() != segid_ || atom.get_iseg() != iseg_) {
            throw std::invalid_argument("Atom does not belong to this residue");
        }
        auto atom_ptr = std::make_shared<Atom>(atom);
        atoms_.push_back(atom_ptr);
        update_atom_map();
    }

    void add_atom(std::shared_ptr<Atom> atom) {
        if (!atom) {
            throw std::invalid_argument("Null atom pointer");
        }
        if (atom->get_resname() != resname_ || atom->get_ires() != ires_ ||
            atom->get_segid() != segid_ || atom->get_iseg() != iseg_) {
            throw std::invalid_argument("Atom does not belong to this residue");
        }
        atoms_.push_back(atom);
        update_atom_map();
    }

    void remove_atom(const std::string& atom_type) {
        auto it = std::remove_if(atoms_.begin(), atoms_.end(),
            [&atom_type](const std::shared_ptr<Atom>& atom) {
                return atom && atom->get_type() == atom_type;
            });
        if (it != atoms_.end()) {
            atoms_.erase(it, atoms_.end());
            update_atom_map();
        }
    }

    void clear_atoms() {
        atoms_.clear();
        atom_map_.clear();
    }

    // Atom access
    const std::vector<std::shared_ptr<Atom>>& get_atoms() const noexcept { 
        return atoms_; 
    }

    std::shared_ptr<Atom> find_atom(const std::string& type) const {
        auto it = atom_map_.find(type);
        return (it != atom_map_.end()) ? it->second : nullptr;
    }

    std::shared_ptr<Atom> find_atom_by_index(size_t index) const {
        return (index < atoms_.size()) ? atoms_[index] : nullptr;
    }

    size_t atom_count() const noexcept { return atoms_.size(); }

    // Atom selection
    std::vector<std::shared_ptr<Atom>> select_atoms(
        const std::function<bool(const Atom&)>& predicate) const {
        std::vector<std::shared_ptr<Atom>> selected;
        for (const auto& atom : atoms_) {
            if (atom && predicate(*atom)) {
                selected.push_back(atom);
            }
        }
        return selected;
    }

    std::vector<std::shared_ptr<Atom>> select_atoms_by_type(
        const std::vector<std::string>& types) const {
        std::vector<std::shared_ptr<Atom>> selected;
        for (const auto& type : types) {
            auto atom = find_atom(type);
            if (atom) selected.push_back(atom);
        }
        return selected;
    }

    bool has_atom_type(const std::string& type) const {
        return atom_map_.find(type) != atom_map_.end();
    }

    // Center of mass calculation
    void calculate_center_of_mass() {
        com_ = {0.0, 0.0, 0.0};
        double total_mass = 0.0;

        for (const auto& atom : atoms_) {
            if (!atom) continue;
            double mass = atom->get_mass();
            const auto& coor = atom->get_coor();
            for (int i = 0; i < 3; ++i) {
                com_[i] += mass * coor[i];
            }
            total_mass += mass;
        }

        if (total_mass > 0.0) {
            for (double& x : com_) x /= total_mass;
        }
    }

    const std::array<double, 3>& get_center_of_mass() const noexcept {
        return com_;
    }

    // Total properties calculation
    double get_total_mass() const {
        double total = 0.0;
        for (const auto& atom : atoms_) {
            if (atom) total += atom->get_mass();
        }
        return total;
    }

    double get_total_charge() const {
        double total = 0.0;
        for (const auto& atom : atoms_) {
            if (atom) total += atom->get_charge();
        }
        return total;
    }

    // IValidatable interface
    bool is_valid() const override {
        if (resname_.empty() || ires_ <= 0 || segid_.empty()) {
            return false;
        }
        
        return std::all_of(atoms_.begin(), atoms_.end(),
                          [](const std::shared_ptr<Atom>& atom) {
                              return atom && atom->is_valid();
                          });
    }

    std::string get_validation_error() const override {
        if (resname_.empty()) return "Empty residue name";
        if (ires_ <= 0) return "Invalid residue number";
        if (segid_.empty()) return "Empty segment ID";
        
        for (const auto& atom : atoms_) {
            if (!atom) return "Null atom pointer";
            if (!atom->is_valid()) return "Invalid atom: " + atom->get_validation_error();
        }
        
        return "";
    }

    // IIdentifiable interface
    std::string get_id() const override {
        return utils::id::generate_residue_id(resname_, ires_, segid_, chain_);
    }

    void set_id(const std::string& id) override {
        // Parse ID format: segid[:chain]:resnum:resname
        std::vector<std::string> parts;
        std::istringstream iss(id);
        std::string part;
        while (std::getline(iss, part, ':')) {
            parts.push_back(part);
        }
        
        if (parts.size() >= 3) {
            size_t idx = 0;
            segid_ = parts[idx++];
            
            // Check if chain is present
            if (parts.size() == 4) {
                chain_ = parts[idx].empty() ? ' ' : parts[idx][0];
                idx++;
            }
            
            ires_ = std::stoi(parts[idx++]);
            resname_ = parts[idx];
        }
    }

protected:
    void update_atom_map() {
        atom_map_.clear();
        for (const auto& atom : atoms_) {
            if (atom) {
                atom_map_[atom->get_type()] = atom;
                // Also store PDB-formatted name
                atom_map_[atom->get_formatted_atom_name()] = atom;
            }
        }
    }

private:
    // CHARMM standard fields
    std::string resname_;                ///< Residue name (RESNAME)
    int ires_;                          ///< Residue number (IRES)
    std::string segid_;                 ///< Segment ID (SEGID)
    int iseg_;                          ///< Segment number (ISEG)
    int igro_;                          ///< Group number (IGRO)
    char chain_;                        ///< Chain identifier
    char inscode_;                      ///< Insertion code
    int move_;                          ///< Movement flag (MOVE)
    int ignore_;                        ///< Ignore flag (IGNORE)
    int constrain_;                     ///< Constraint flag (CONSTRAIN)
    bool hetatm_;                       ///< HETATM flag

    // Atom storage
    std::vector<std::shared_ptr<Atom>> atoms_;      ///< Atoms in residue
    std::unordered_map<std::string, std::shared_ptr<Atom>> atom_map_;  ///< Quick atom lookup

    // Secondary structure
    SecondaryStructure sec_struct_;     ///< Secondary structure type

    // Calculated properties
    std::array<double, 3> com_{0.0, 0.0, 0.0};  ///< Center of mass

    // Scalar properties (SCA1-9)
    std::array<double, 9> scalar_;      ///< User-defined scalar properties
};

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_RESIDUE_COMPOSITE_HPP 