#pragma once

#ifndef PYGCMC_MODEL_MOLECULAR_COMPOSITE_HPP
#define PYGCMC_MODEL_MOLECULAR_COMPOSITE_HPP

#include <vector>
#include <memory>
#include <map>
#include <string>
#include <unordered_map>
#include <array>
#include <set>
#include <tuple>
#include <algorithm>
#include "../atom/AtomMain.hpp"
#include "../residue/ResidueMain.hpp"
#include "../common/ModelInterface.hpp"
#include "../common/ModelUtils.hpp"

namespace pygcmc {
namespace model {

/**
 * @brief Standardized CMAP structure for unified PSF and TOP format processing
 */
struct StandardCmap {
    std::array<int, 5> atoms;     ///< Standardized 5 atom indices
    std::array<int, 8> raw_atoms; ///< Original format atom indices (8 for PSF, 5+3 of -1 for TOP)
    bool is_psf_format;           ///< Whether it is PSF format
    int function_type = 1;        ///< CMAP function type
    
    StandardCmap() : is_psf_format(false) {
        atoms.fill(-1);
        raw_atoms.fill(-1);
    }
};

/**
 * @brief Terminal information structure
 */
struct TerminalInfo {
    std::string chain_id;
    int start_residue = 0;
    int end_residue = 0;
    std::string terminal_type;  // "N-terminal", "C-terminal", etc.
};

/**
 * @brief Secondary structure information
 */
struct SecondaryStructure {
    std::string type;           // "HELIX", "SHEET", etc.
    std::string identifier;
    int start_residue = 0;
    int end_residue = 0;
    std::string chain_id;
    std::string comment;
};

/**
 * @brief Core molecular functionality for residue composition and topology management
 */
class MolecularComposite : public IValidatable, public IIdentifiable {
public:
    // Constructors
    MolecularComposite() = default;
    
    explicit MolecularComposite(const std::string& name) : name_(name) {}

    // Basic properties
    const std::string& get_name() const noexcept { return name_; }
    void set_name(const std::string& name) { name_ = name; }

    // Atom management
    void add_atom(std::shared_ptr<Atom> atom) {
        if (!atom) {
            throw std::invalid_argument("Cannot add null atom");
        }
        atoms_.push_back(atom);
        update_atom_map();
    }

    void remove_atom(size_t index) {
        if (index < atoms_.size()) {
            atoms_.erase(atoms_.begin() + index);
            update_atom_map();
        }
    }

    const std::vector<std::shared_ptr<Atom>>& get_atoms() const noexcept {
        return atoms_;
    }

    std::shared_ptr<Atom> get_atom(size_t index) const {
        return (index < atoms_.size()) ? atoms_[index] : nullptr;
    }

    size_t get_num_atoms() const noexcept { return atoms_.size(); }

    // Residue management
    void add_residue(std::shared_ptr<Residue> residue) {
        if (!residue) {
            throw std::invalid_argument("Cannot add null residue");
        }
        residues_.push_back(residue);
        update_residue_map();
    }

    void remove_residue(size_t index) {
        if (index < residues_.size()) {
            residues_.erase(residues_.begin() + index);
            update_residue_map();
        }
    }

    const std::vector<std::shared_ptr<Residue>>& get_residues() const noexcept {
        return residues_;
    }

    std::shared_ptr<Residue> get_residue(size_t index) const {
        return (index < residues_.size()) ? residues_[index] : nullptr;
    }

    size_t get_num_residues() const noexcept { return residues_.size(); }

    // Atom lookup methods
    std::shared_ptr<Atom> find_atom_by_id(int atom_id) const {
        auto it = std::find_if(atoms_.begin(), atoms_.end(),
            [atom_id](const std::shared_ptr<Atom>& atom) {
                return atom && atom->get_bynu() == atom_id;
            });
        return (it != atoms_.end()) ? *it : nullptr;
    }

    std::shared_ptr<Atom> find_atom_by_name(const std::string& residue_name, 
                                           int residue_number, 
                                           const std::string& atom_name) const {
        auto key = std::make_tuple(residue_name, residue_number, atom_name);
        auto it = atom_map_.find(key);
        return (it != atom_map_.end() && it->second < atoms_.size()) 
                ? atoms_[it->second] : nullptr;
    }

    // Residue lookup methods
    std::shared_ptr<Residue> find_residue_by_id(const std::string& residue_name, 
                                               int residue_number) const {
        auto key = std::make_pair(residue_name, residue_number);
        auto it = residue_map_.find(key);
        return (it != residue_map_.end() && it->second < residues_.size()) 
                ? residues_[it->second] : nullptr;
    }

    std::vector<std::shared_ptr<Residue>> find_residues_by_name(const std::string& name) const {
        std::vector<std::shared_ptr<Residue>> result;
        for (const auto& residue : residues_) {
            if (residue && residue->get_resname() == name) {
                result.push_back(residue);
            }
        }
        return result;
    }

    // Chain management
    void add_chain(const std::string& chain_id) {
        if (std::find(chain_ids_.begin(), chain_ids_.end(), chain_id) == chain_ids_.end()) {
            chain_ids_.push_back(chain_id);
        }
    }

    const std::vector<std::string>& get_chain_ids() const noexcept {
        return chain_ids_;
    }

    std::vector<std::shared_ptr<Residue>> get_residues_by_chain(char chain_id) const {
        std::vector<std::shared_ptr<Residue>> result;
        for (const auto& residue : residues_) {
            if (residue && residue->get_chain() == chain_id) {
                result.push_back(residue);
            }
        }
        return result;
    }

    // Box dimensions
    void set_box_dimensions(const std::vector<double>& dimensions) {
        box_dimensions_ = dimensions;
    }

    const std::vector<double>& get_box_dimensions() const noexcept {
        return box_dimensions_;
    }

    // Terminal information
    void add_terminal(const TerminalInfo& terminal) {
        terminals_.push_back(terminal);
    }

    const std::vector<TerminalInfo>& get_terminals() const noexcept {
        return terminals_;
    }

    // Secondary structure
    void add_secondary_structure(const std::string& type, const SecondaryStructure& ss) {
        secondary_structures_[type].push_back(ss);
    }

    const std::map<std::string, std::vector<SecondaryStructure>>& 
    get_secondary_structures() const noexcept {
        return secondary_structures_;
    }

    // CMAP management
    void add_cmap(const StandardCmap& cmap) {
        cmaps_.push_back(cmap);
    }

    const std::vector<StandardCmap>& get_cmaps() const noexcept {
        return cmaps_;
    }

    size_t get_num_cmaps() const noexcept { return cmaps_.size(); }

    // Title information
    void add_title(const std::string& title) {
        titles_.push_back(title);
    }

    const std::vector<std::string>& get_titles() const noexcept {
        return titles_;
    }

    // Disulfide bonds
    void add_disulfide_bond(const std::string& bond) {
        disulfide_bonds_.push_back(bond);
    }

    const std::vector<std::string>& get_disulfide_bonds() const noexcept {
        return disulfide_bonds_;
    }

    // Calculate molecular properties
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

    std::array<double, 3> get_center_of_mass() const {
        std::array<double, 3> com = {0.0, 0.0, 0.0};
        double total_mass = 0.0;

        for (const auto& atom : atoms_) {
            if (!atom) continue;
            double mass = atom->get_mass();
            const auto& coor = atom->get_coor();
            for (int i = 0; i < 3; ++i) {
                com[i] += mass * coor[i];
            }
            total_mass += mass;
        }

        if (total_mass > 0.0) {
            for (double& x : com) x /= total_mass;
        }

        return com;
    }

    // Clear all data
    void clear() {
        atoms_.clear();
        residues_.clear();
        chain_ids_.clear();
        terminals_.clear();
        secondary_structures_.clear();
        disulfide_bonds_.clear();
        cmaps_.clear();
        titles_.clear();
        box_dimensions_.clear();
        
        atom_map_.clear();
        residue_map_.clear();
    }

    // IValidatable interface
    bool is_valid() const override {
        // Check all atoms are valid
        for (const auto& atom : atoms_) {
            if (!atom || !atom->is_valid()) {
                return false;
            }
        }
        
        // Check all residues are valid
        for (const auto& residue : residues_) {
            if (!residue || !residue->is_valid()) {
                return false;
            }
        }
        
        return true;
    }

    std::string get_validation_error() const override {
        for (size_t i = 0; i < atoms_.size(); ++i) {
            if (!atoms_[i]) {
                return "Null atom at index " + std::to_string(i);
            }
            if (!atoms_[i]->is_valid()) {
                return "Invalid atom at index " + std::to_string(i) + ": " + 
                       atoms_[i]->get_validation_error();
            }
        }
        
        for (size_t i = 0; i < residues_.size(); ++i) {
            if (!residues_[i]) {
                return "Null residue at index " + std::to_string(i);
            }
            if (!residues_[i]->is_valid()) {
                return "Invalid residue at index " + std::to_string(i) + ": " + 
                       residues_[i]->get_validation_error();
            }
        }
        
        return "";
    }

    // IIdentifiable interface
    std::string get_id() const override {
        return name_.empty() ? "molecular_" + std::to_string(atoms_.size()) : name_;
    }

    void set_id(const std::string& id) override {
        name_ = id;
    }

protected:
    void update_atom_map() {
        atom_map_.clear();
        for (size_t i = 0; i < atoms_.size(); ++i) {
            const auto& atom = atoms_[i];
            if (atom) {
                auto key = std::make_tuple(atom->get_resname(), 
                                         atom->get_ires(), 
                                         atom->get_type());
                atom_map_[key] = i;
            }
        }
    }

    void update_residue_map() {
        residue_map_.clear();
        for (size_t i = 0; i < residues_.size(); ++i) {
            const auto& residue = residues_[i];
            if (residue) {
                auto key = std::make_pair(residue->get_resname(), residue->get_ires());
                residue_map_[key] = i;
            }
        }
    }

private:
    std::string name_;                                      ///< Molecular name/identifier

    // Structure data
    std::vector<std::shared_ptr<Atom>> atoms_;              ///< Atom list
    std::vector<std::shared_ptr<Residue>> residues_;        ///< Residue list
    std::vector<std::string> chain_ids_;                    ///< Chain identifiers
    std::vector<TerminalInfo> terminals_;                   ///< Terminal information
    std::map<std::string, std::vector<SecondaryStructure>> secondary_structures_;  ///< Secondary structures
    std::vector<std::string> disulfide_bonds_;              ///< Disulfide bonds
    std::vector<double> box_dimensions_;                    ///< Box dimensions

    // CMAP data
    std::vector<StandardCmap> cmaps_;                       ///< CMAP information
    
    // Title information
    std::vector<std::string> titles_;                       ///< Title information

    // Lookup maps
    std::map<std::tuple<std::string, int, std::string>, size_t> atom_map_;     ///< (resname, resnum, atomname) -> index
    std::map<std::pair<std::string, int>, size_t> residue_map_;               ///< (resname, resnum) -> index
};

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MOLECULAR_COMPOSITE_HPP 