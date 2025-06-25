#pragma once

#ifndef PYGCMC_MODEL_MOLECULE_MAIN_HPP
#define PYGCMC_MODEL_MOLECULE_MAIN_HPP

#include <vector>
#include <memory>
#include <map>
#include <string>
#include <unordered_map>
#include <array>
#include <set>
#include "../atom/AtomMain.hpp"
#include "../residue/ResidueMain.hpp"
#include "../structure/StructureMain.hpp"
#include "../topology/TopologyMain.hpp"

namespace pygcmc {
namespace model {
namespace molecule {

/**
 * @brief Standardized CMAP structure, used for unified processing of PSF and TOP formats
 */
struct StandardCmap {
    std::array<int, 5> atoms;     ///< Standardized 5 atom indices
    std::array<int, 8> raw_atoms; ///< Original format atom indices (8 for PSF, 5+3 of -1 for TOP)
    bool is_psf_format;           ///< Whether it is PSF format
    int function_type = 1;        ///< CMAP function type
};

/**
 * @brief Molecular class, merges Structure and Topology data for subsequent calculations
 */
class Molecular {
public:
    Molecular() = default;
    ~Molecular() = default;

    // Structure (PDB) related data
    std::vector<std::shared_ptr<atom::Atom>> atoms;  // Atom list
    std::vector<std::shared_ptr<residue::Residue>> residues;  // Residue list
    std::vector<structure::Structure::TerminalInfo> terminals;  // Chain termination information
    std::map<std::string, std::vector<structure::Structure::SecondaryStructure>> helices;  // Helical structures
    std::map<std::string, std::vector<std::string>> sheets;  // Beta-sheet structures
    std::vector<std::string> ssbonds;  // Disulfide bonds
    std::vector<double> boxDimensions;  // Box dimensions

    // Topology (PSF/TOP) related data
    std::vector<topology::TopologyAtom> topology_atoms;  // Atom information in topology (includes charge, mass, etc.)
    std::vector<topology::TopologyResidue> topology_residues;  // Residue information in topology
    std::vector<topology::TopologySegment> segments;  // Fragment information
    std::vector<topology::TopologyBond> bonds;  // Bond
    std::vector<topology::TopologyAngle> angles;  // Bond angle
    std::vector<topology::TopologyDihedral> dihedrals;  // Dihedral angle (including improper)
    std::vector<topology::TopologyDonor> donors;  // Hydrogen bond donor
    std::vector<topology::TopologyAcceptor> acceptors;  // Hydrogen bond acceptor
    std::map<int, std::set<int>> exclusions;  // Non-bond exclusion
    std::vector<topology::TopologyGroup> groups;  // Atom group
    std::vector<topology::TopologyCmap> cmaps;  // Original CMAP item
    std::vector<StandardCmap> standard_cmaps;  // Standardized CMAP item
    std::vector<std::string> titles;  // PSF file title information

    // Lookup mapping
    std::unordered_map<std::string, int> segment_map;  // segment_name -> index
    std::map<std::pair<std::string, int>, int> residue_map;  // (residue_name, number) -> index
    std::map<std::tuple<std::string, int, std::string>, int> atom_map;  // (residue_name, number, atom_name) -> index

    // Get atom count
    size_t get_num_atoms() const { return atoms.size(); }
    
    // Get residue count
    size_t get_num_residues() const { return residues.size(); }
    
    // Get fragment count
    size_t get_num_segments() const { return segments.size(); }
    
    // Get bond count
    size_t get_num_bonds() const { return bonds.size(); }
    
    // Get bond angle count
    size_t get_num_angles() const { return angles.size(); }
    
    // Get dihedral angle count (excluding improper)
    size_t get_num_dihedrals() const {
        size_t count = 0;
        for (const auto& dihedral : dihedrals) {
            if (!dihedral.improper) count++;
        }
        return count;
    }
    
    // Get improper count
    size_t get_num_impropers() const {
        size_t count = 0;
        for (const auto& dihedral : dihedrals) {
            if (dihedral.improper) count++;
        }
        return count;
    }

    // Standardized CMAP related methods
    void add_standard_cmap(const topology::TopologyCmap& cmap) {
        StandardCmap std_cmap;
        std_cmap.raw_atoms = cmap.atoms;
        std_cmap.is_psf_format = (cmap.atoms[5] != -1);  // Determine if it is PSF format
        
        // Set standardized 5 atoms
        if (std_cmap.is_psf_format) {
            // PSF format: Use the first 4 atoms and the 8th atom
            for (int i = 0; i < 4; ++i) {
                std_cmap.atoms[i] = cmap.atoms[i];
            }
            std_cmap.atoms[4] = cmap.atoms[7];  // Use the 8th atom as the 5th atom
        } else {
            // TOP format: Directly use the first 5 atoms
            for (int i = 0; i < 5; ++i) {
                std_cmap.atoms[i] = cmap.atoms[i];
            }
        }
        std_cmap.function_type = cmap.function_type;
        standard_cmaps.push_back(std_cmap);
    }

    // Get standardized CMAP count
    size_t get_num_standard_cmaps() const { return standard_cmaps.size(); }

    // Clear all data
    void clear() {
        // Structure data
        atoms.clear();
        residues.clear();
        terminals.clear();
        helices.clear();
        sheets.clear();
        ssbonds.clear();
        boxDimensions.clear();

        // Topology data
        topology_atoms.clear();
        topology_residues.clear();
        segments.clear();
        bonds.clear();
        angles.clear();
        dihedrals.clear();
        donors.clear();
        acceptors.clear();
        exclusions.clear();
        groups.clear();
        cmaps.clear();
        standard_cmaps.clear();
        titles.clear();

        // Lookup mapping
        segment_map.clear();
        residue_map.clear();
        atom_map.clear();
    }
};

} // namespace molecule

// Backward compatibility: provide Molecular in the model namespace
using Molecular = molecule::Molecular;

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_MOLECULE_MAIN_HPP