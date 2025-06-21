#pragma once

#ifndef PYGCMC_MODEL_MOLECULE_COMPOSITE_HPP
#define PYGCMC_MODEL_MOLECULE_COMPOSITE_HPP

#include "../common/ModelInterface.hpp"
#include "../atom/AtomMain.hpp"
#include "../residue/ResidueMain.hpp"
#include "../topology/TopologyMain.hpp"
#include "../structure.hpp"
#include <vector>
#include <memory>
#include <map>
#include <string>
#include <unordered_map>
#include <array>
#include <set>
#include <tuple>

namespace pygcmc {
namespace model {
namespace molecule {

/**
 * @brief Standardized CMAP structure for unified PSF and TOP format processing
 */
struct StandardCmap {
    std::array<int, 5> atoms;     ///< Standardized 5 atom indices
    std::array<int, 8> raw_atoms; ///< Original format atom indices (8 for PSF, 5+3 of -1 for TOP)
    bool is_psf_format;           ///< Whether it is PSF format
    int function_type = 1;        ///< CMAP function type
    std::string description;      ///< Description information
};

/**
 * @brief Structure information aggregator
 */
struct StructureInfo {
    // Use original Structure types directly for full compatibility
    using TerminalInfo = pygcmc::model::Structure::TerminalInfo;
    using SecondaryStructure = pygcmc::model::Structure::SecondaryStructure;
    
    std::vector<TerminalInfo> terminals;
    std::map<std::string, std::vector<SecondaryStructure>> helices;
    std::map<std::string, std::vector<std::string>> sheets;
    std::vector<std::string> ssbonds;
    std::vector<double> box_dimensions;
};

/**
 * @brief Topology information from PSF/TOP files - using original topology types
 */
struct TopologyInfo {
    // Use original topology types for compatibility
    std::vector<topology::TopologyAtom> atoms;
    std::vector<topology::TopologyResidue> residues;
    std::vector<topology::TopologySegment> segments;
    std::vector<topology::TopologyBond> bonds;
    std::vector<topology::TopologyAngle> angles;
    std::vector<topology::TopologyDihedral> dihedrals;
    std::vector<topology::TopologyDonor> donors;
    std::vector<topology::TopologyAcceptor> acceptors;
    std::map<int, std::set<int>> exclusions;
    std::vector<topology::TopologyGroup> groups;
    std::vector<topology::TopologyCmap> cmaps;
    std::vector<std::string> titles;
};

/**
 * @brief Core molecular composite class that manages atoms and residues
 * This class handles the composition of molecular systems from structure and topology data
 */
class MolecularComposite : public common::IValidatable {
public:
    MolecularComposite() = default;
    virtual ~MolecularComposite() = default;

    // IValidatable interface
    bool is_valid() const override {
        // Check basic consistency
        if (atoms.empty()) return false;
        
        // Check all atoms are valid
        for (const auto& atom : atoms) {
            if (!atom || !atom->is_valid()) return false;
        }
        
        // Check all residues are valid
        for (const auto& residue : residues) {
            if (!residue || !residue->is_valid()) return false;
        }
        
        return true;
    }

    // Core data access
    const std::vector<std::shared_ptr<atom::Atom>>& get_atoms() const { return atoms; }
    const std::vector<std::shared_ptr<residue::Residue>>& get_residues() const { return residues; }
    const StructureInfo& get_structure_info() const { return structure_info; }
    const TopologyInfo& get_topology_info() const { return topology_info; }

    std::vector<std::shared_ptr<atom::Atom>>& get_atoms() { return atoms; }
    std::vector<std::shared_ptr<residue::Residue>>& get_residues() { return residues; }
    StructureInfo& get_structure_info() { return structure_info; }
    TopologyInfo& get_topology_info() { return topology_info; }

    // Size getters
    size_t get_num_atoms() const { return atoms.size(); }
    size_t get_num_residues() const { return residues.size(); }
    size_t get_num_segments() const { return topology_info.segments.size(); }
    size_t get_num_bonds() const { return topology_info.bonds.size(); }
    size_t get_num_angles() const { return topology_info.angles.size(); }

    size_t get_num_dihedrals() const {
        size_t count = 0;
        for (const auto& dihedral : topology_info.dihedrals) {
            if (!dihedral.improper) count++;
        }
        return count;
    }

    size_t get_num_impropers() const {
        size_t count = 0;
        for (const auto& dihedral : topology_info.dihedrals) {
            if (dihedral.improper) count++;
        }
        return count;
    }

    size_t get_num_standard_cmaps() const { return standard_cmaps.size(); }

    // Atom management
    void add_atom(std::shared_ptr<atom::Atom> atom) {
        if (!atom) throw std::invalid_argument("Null atom pointer");
        atoms.push_back(atom);
        refresh_atom_map();
    }

    void add_atoms(const std::vector<std::shared_ptr<atom::Atom>>& new_atoms) {
        for (const auto& atom : new_atoms) {
            if (atom) atoms.push_back(atom);
        }
        refresh_atom_map();
    }

    void remove_atom(size_t index) {
        if (index < atoms.size()) {
            atoms.erase(atoms.begin() + index);
            refresh_atom_map();
        }
    }

    std::shared_ptr<atom::Atom> find_atom(const std::string& segment_id, 
                                          [[maybe_unused]] const std::string& residue_name, 
                                          int residue_number,
                                          const std::string& atom_name) const {
        auto key = std::make_tuple(segment_id, residue_number, atom_name);
        auto it = atom_map.find(key);
        return (it != atom_map.end() && it->second < static_cast<int>(atoms.size())) ? 
               atoms[it->second] : nullptr;
    }

    // Residue management
    void add_residue(std::shared_ptr<residue::Residue> residue) {
        if (!residue) throw std::invalid_argument("Null residue pointer");
        residues.push_back(residue);
        update_residue_map();
    }

    void add_residues(const std::vector<std::shared_ptr<residue::Residue>>& new_residues) {
        for (const auto& residue : new_residues) {
            if (residue) residues.push_back(residue);
        }
        update_residue_map();
    }

    void remove_residue(size_t index) {
        if (index < residues.size()) {
            residues.erase(residues.begin() + index);
            update_residue_map();
        }
    }

    std::shared_ptr<residue::Residue> find_residue(const std::string& residue_name, 
                                                   int residue_number) const {
        auto key = std::make_pair(residue_name, residue_number);
        auto it = residue_map.find(key);
        return (it != residue_map.end() && it->second < static_cast<int>(residues.size())) ? 
               residues[it->second] : nullptr;
    }

    // CMAP management  
    void add_cmap(const topology::TopologyCmap& cmap) {
        topology_info.cmaps.push_back(cmap);
        add_standard_cmap(cmap);
    }

    void add_standard_cmap(const topology::TopologyCmap& cmap) {
        StandardCmap std_cmap;
        std_cmap.raw_atoms = cmap.atoms;
        std_cmap.is_psf_format = (cmap.atoms[5] != -1);
        
        // Set standardized 5 atoms
        if (std_cmap.is_psf_format) {
            // PSF format: Use the first 4 atoms and the 8th atom
            for (int i = 0; i < 4; ++i) {
                std_cmap.atoms[i] = cmap.atoms[i];
            }
            std_cmap.atoms[4] = cmap.atoms[7];
        } else {
            // TOP format: Directly use the first 5 atoms
            for (int i = 0; i < 5; ++i) {
                std_cmap.atoms[i] = cmap.atoms[i];
            }
        }
        std_cmap.function_type = cmap.function_type;
        standard_cmaps.push_back(std_cmap);
    }

    const std::vector<StandardCmap>& get_standard_cmaps() const { return standard_cmaps; }
    std::vector<StandardCmap>& get_standard_cmaps() { return standard_cmaps; }

    // Clear all data
    void clear() {
        // Structure data
        atoms.clear();
        residues.clear();
        structure_info = StructureInfo{};
        
        // Topology data
        topology_info = TopologyInfo{};
        standard_cmaps.clear();

        // Lookup mappings
        segment_map.clear();
        residue_map.clear();
        atom_map.clear();
    }

    // Mapping utilities
    const std::unordered_map<std::string, int>& get_segment_map() const { return segment_map; }
    const std::map<std::pair<std::string, int>, int>& get_residue_map() const { return residue_map; }
    const std::map<std::tuple<std::string, int, std::string>, int>& get_atom_map() const { return atom_map; }

protected:
    void refresh_atom_map() {
        atom_map.clear();
        for (size_t i = 0; i < atoms.size(); ++i) {
            const auto& atom = atoms[i];
            if (atom) {
                auto key = std::make_tuple(atom->get_segid(), atom->get_ires(), atom->get_type());
                atom_map[key] = static_cast<int>(i);
            }
        }
    }

    void update_residue_map() {
        residue_map.clear();
        for (size_t i = 0; i < residues.size(); ++i) {
            const auto& residue = residues[i];
            if (residue) {
                auto key = std::make_pair(residue->get_resname(), residue->get_ires());
                residue_map[key] = static_cast<int>(i);
            }
        }
    }

    void update_segment_map() {
        segment_map.clear();
        for (size_t i = 0; i < topology_info.segments.size(); ++i) {
            const auto& segment = topology_info.segments[i];
            segment_map[segment.name] = static_cast<int>(i);
        }
    }

private:
    // Core molecular data
    std::vector<std::shared_ptr<atom::Atom>> atoms;
    std::vector<std::shared_ptr<residue::Residue>> residues;
    
    // Structure and topology information
    StructureInfo structure_info;
    TopologyInfo topology_info;
    
    // Standardized CMAP data
    std::vector<StandardCmap> standard_cmaps;

    // Lookup mappings for fast access
    std::unordered_map<std::string, int> segment_map;  // segment_name -> index
    std::map<std::pair<std::string, int>, int> residue_map;  // (residue_name, number) -> index
    std::map<std::tuple<std::string, int, std::string>, int> atom_map;  // (segid, resnum, atom_name) -> index
};

} // namespace molecule
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MOLECULE_COMPOSITE_HPP 