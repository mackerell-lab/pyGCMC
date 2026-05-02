#pragma once

#include "TopologyStructures.hpp"
#include <algorithm>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Query operations for topology elements
 */
class TopologyQueries {
public:
    // Find operations
    static inline std::optional<int> findAtom(
        const std::vector<TopologySegment>& segments,
        const std::map<std::tuple<std::string, int, std::string, std::string>, int>& atom_map,
        const std::string& residue_name, int residue_number, const std::string& atom_name
    ) {
        for (const auto& segment : segments) {
            auto it = atom_map.find(std::make_tuple(residue_name, residue_number, segment.name, atom_name));
            if (it != atom_map.end()) {
                return it->second;
            }
        }
        return std::nullopt;
    }

    static inline std::optional<int> findResidue(
        const std::vector<TopologySegment>& segments,
        const std::map<std::tuple<std::string, int, std::string>, int>& residue_map,
        const std::string& name, int number
    ) {
        for (const auto& segment : segments) {
            auto it = residue_map.find(std::make_tuple(name, number, segment.name));
            if (it != residue_map.end()) {
                return it->second;
            }
        }
        return std::nullopt;
    }

    static inline std::optional<int> findSegment(
        const std::unordered_map<std::string, int>& segment_map,
        const std::string& name
    ) {
        auto it = segment_map.find(name);
        return (it != segment_map.end()) ? std::optional<int>(it->second) : std::nullopt;
    }

    // Check connectivity operations
    static inline bool hasBond(const std::vector<TopologyBond>& bonds, int atom1, int atom2) {
        for (const auto& bond : bonds) {
            if ((bond.atom1 == atom1 && bond.atom2 == atom2) ||
                (bond.atom1 == atom2 && bond.atom2 == atom1)) {
                return true;
            }
        }
        return false;
    }

    static inline bool hasAngle(const std::vector<TopologyAngle>& angles, int atom1, int atom2, int atom3) {
        for (const auto& angle : angles) {
            if ((angle.atom1 == atom1 && angle.atom2 == atom2 && angle.atom3 == atom3) ||
                (angle.atom1 == atom3 && angle.atom2 == atom2 && angle.atom3 == atom1)) {
                return true;
            }
        }
        return false;
    }

    static inline bool hasDihedral(const std::vector<TopologyDihedral>& dihedrals,
                                  int atom1, int atom2, int atom3, int atom4) {
        for (const auto& dihedral : dihedrals) {
            if (dihedral.improper) continue;
            if ((dihedral.atom1 == atom1 && dihedral.atom2 == atom2 &&
                 dihedral.atom3 == atom3 && dihedral.atom4 == atom4) ||
                (dihedral.atom1 == atom4 && dihedral.atom2 == atom3 &&
                 dihedral.atom3 == atom2 && dihedral.atom4 == atom1)) {
                return true;
            }
        }
        return false;
    }

    static inline bool hasImproper(const std::vector<TopologyDihedral>& dihedrals,
                                  int atom1, int atom2, int atom3, int atom4) {
        for (const auto& dihedral : dihedrals) {
            if (!dihedral.improper) continue;
            if ((dihedral.atom1 == atom1 && dihedral.atom2 == atom2 &&
                 dihedral.atom3 == atom3 && dihedral.atom4 == atom4) ||
                (dihedral.atom1 == atom4 && dihedral.atom2 == atom3 &&
                 dihedral.atom3 == atom2 && dihedral.atom4 == atom1)) {
                return true;
            }
        }
        return false;
    }

    // Check special features
    static inline bool hasDonor(const std::vector<TopologyDonor>& donors, int donor_atom, int hydrogen_atom) {
        for (const auto& donor : donors) {
            if (donor.donor_atom == donor_atom && donor.hydrogen_atom == hydrogen_atom) {
                return true;
            }
        }
        return false;
    }

    static inline bool hasAcceptor(const std::vector<TopologyAcceptor>& acceptors, int acceptor_atom) {
        for (const auto& acceptor : acceptors) {
            if (acceptor.acceptor_atom == acceptor_atom) {
                return true;
            }
        }
        return false;
    }

    static inline bool hasGroup(const std::vector<TopologyGroup>& groups, int group_id) {
        for (const auto& group : groups) {
            if (group.id == group_id) {
                return true;
            }
        }
        return false;
    }

    static inline bool hasCmap(const std::vector<TopologyCmap>& cmaps) {
        return !cmaps.empty();
    }

    static inline bool hasCmap(const std::vector<TopologyCmap>& cmaps, const std::array<int, 8>& atoms) {
        for (const auto& cmap : cmaps) {
            if (cmap.atoms == atoms) return true;
        }
        return false;
    }

    static inline bool hasCmap(const std::vector<TopologyCmap>& cmaps, const std::vector<int>& atoms) {
        if (atoms.size() == 8) {
            for (const auto& cmap : cmaps) {
                bool match = true;
                for (size_t i = 0; i < 8; ++i) {
                    if (cmap.atoms[i] != atoms[i]) {
                        match = false;
                        break;
                    }
                }
                if (match) return true;
            }
        }
        else if (atoms.size() == 5) {
            for (const auto& cmap : cmaps) {
                bool match = true;
                for (size_t i = 0; i < 5; ++i) {
                    if (cmap.atoms[i] != atoms[i]) {
                        match = false;
                        break;
                    }
                }
                if (match) return true;
            }
        }
        return false;
    }

    // Count operations
    static inline size_t countDihedrals(const std::vector<TopologyDihedral>& dihedrals) {
        return std::count_if(dihedrals.begin(), dihedrals.end(),
                            [](const TopologyDihedral& d) { return !d.improper; });
    }

    static inline size_t countImpropers(const std::vector<TopologyDihedral>& dihedrals) {
        return std::count_if(dihedrals.begin(), dihedrals.end(),
                            [](const TopologyDihedral& d) { return d.improper; });
    }

    // Utility checks
    static inline bool hasAtom(const std::vector<TopologyAtom>& atoms, int index) {
        return index >= 0 && index < static_cast<int>(atoms.size());
    }

    static inline bool hasResidue(const std::vector<TopologyResidue>& residues, int index) {
        return index >= 0 && index < static_cast<int>(residues.size());
    }

    static inline bool hasSegment(const std::vector<TopologySegment>& segments, int index) {
        return index >= 0 && index < static_cast<int>(segments.size());
    }
};

} // namespace topology
} // namespace model
} // namespace pygcmc
