#pragma once

#ifndef PYGCMC_MODEL_TOPOLOGY_MAIN_IMPL_HPP
#define PYGCMC_MODEL_TOPOLOGY_MAIN_IMPL_HPP

#include "TopologyMain.hpp"
#include <algorithm>

namespace pygcmc {
namespace model {
namespace topology {

// Implementation of find methods
inline std::optional<int> Topology::find_atom(const std::string& residue_name, int residue_number, 
                                              const std::string& atom_name) const {
    for (const auto& segment : storage_.segments) {
        auto it = storage_.atom_map.find(std::make_tuple(residue_name, residue_number, segment.name, atom_name));
        if (it != storage_.atom_map.end()) {
            return it->second;
        }
    }
    return std::nullopt;
}

inline std::optional<int> Topology::find_residue(const std::string& name, int number) const {
    for (const auto& segment : storage_.segments) {
        auto it = storage_.residue_map.find(std::make_tuple(name, number, segment.name));
        if (it != storage_.residue_map.end()) {
            return it->second;
        }
    }
    return std::nullopt;
}

inline std::optional<int> Topology::find_segment(const std::string& name) const {
    auto it = storage_.segment_map.find(name);
    if (it != storage_.segment_map.end()) {
        return it->second;
    }
    return std::nullopt;
}

// Implementation of count methods
inline size_t Topology::get_num_dihedrals() const {
    size_t count = 0;
    for (const auto& dih : storage_.dihedrals) {
        if (!dih.improper) count++;
    }
    return count;
}

inline size_t Topology::get_num_impropers() const {
    size_t count = 0;
    for (const auto& dih : storage_.dihedrals) {
        if (dih.improper) count++;
    }
    return count;
}

// Implementation of check methods
inline bool Topology::has_cmap(const std::vector<int>& atoms) const {
    if (atoms.size() < 5) return false;
    return std::any_of(storage_.cmaps.begin(), storage_.cmaps.end(),
        [&atoms](const TopologyCmap& cmap) {
            for (size_t i = 0; i < 5 && i < atoms.size(); ++i) {
                if (cmap.atoms[i] != atoms[i]) return false;
            }
            return true;
        });
}

inline bool Topology::has_bond(int atom1, int atom2) const {
    return std::any_of(storage_.bonds.begin(), storage_.bonds.end(),
        [atom1, atom2](const TopologyBond& bond) {
            return (bond.atom1 == atom1 && bond.atom2 == atom2) ||
                   (bond.atom1 == atom2 && bond.atom2 == atom1);
        });
}

inline bool Topology::has_angle(int atom1, int atom2, int atom3) const {
    return std::any_of(storage_.angles.begin(), storage_.angles.end(),
        [atom1, atom2, atom3](const TopologyAngle& angle) {
            return (angle.atom1 == atom1 && angle.atom2 == atom2 && angle.atom3 == atom3) ||
                   (angle.atom1 == atom3 && angle.atom2 == atom2 && angle.atom3 == atom1);
        });
}

inline bool Topology::has_dihedral(int atom1, int atom2, int atom3, int atom4) const {
    return std::any_of(storage_.dihedrals.begin(), storage_.dihedrals.end(),
        [atom1, atom2, atom3, atom4](const TopologyDihedral& dihedral) {
            return !dihedral.improper &&
                   ((dihedral.atom1 == atom1 && dihedral.atom2 == atom2 && 
                     dihedral.atom3 == atom3 && dihedral.atom4 == atom4) ||
                    (dihedral.atom1 == atom4 && dihedral.atom2 == atom3 && 
                     dihedral.atom3 == atom2 && dihedral.atom4 == atom1));
        });
}

inline bool Topology::has_improper(int atom1, int atom2, int atom3, int atom4) const {
    return std::any_of(storage_.dihedrals.begin(), storage_.dihedrals.end(),
        [atom1, atom2, atom3, atom4](const TopologyDihedral& dihedral) {
            return dihedral.improper &&
                   dihedral.atom1 == atom1 && dihedral.atom2 == atom2 && 
                   dihedral.atom3 == atom3 && dihedral.atom4 == atom4;
        });
}

inline bool Topology::has_donor(int donor_atom) const {
    return std::any_of(storage_.donors.begin(), storage_.donors.end(),
        [donor_atom](const TopologyDonor& d) { 
            return d.donor_atom == donor_atom; 
        });
}

inline bool Topology::has_donor(int donor_atom, int hydrogen_atom) const {
    return std::any_of(storage_.donors.begin(), storage_.donors.end(),
        [donor_atom, hydrogen_atom](const TopologyDonor& d) { 
            return d.donor_atom == donor_atom && d.hydrogen_atom == hydrogen_atom; 
        });
}

inline bool Topology::has_acceptor(int acceptor_atom) const {
    return std::any_of(storage_.acceptors.begin(), storage_.acceptors.end(),
        [acceptor_atom](const TopologyAcceptor& a) { 
            return a.acceptor_atom == acceptor_atom; 
        });
}

inline bool Topology::has_group(int group_id) const {
    return std::any_of(storage_.groups.begin(), storage_.groups.end(),
        [group_id](const TopologyGroup& g) { 
            return g.id == group_id; 
        });
}

} // namespace topology
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_TOPOLOGY_MAIN_IMPL_HPP 