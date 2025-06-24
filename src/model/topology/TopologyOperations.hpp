#pragma once

#include "TopologyStructures.hpp"
#include <stdexcept>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Operations for adding topology elements
 */
class TopologyOperations {
public:
    // Add atom operation
    static inline int addAtom(
        std::vector<TopologyAtom>& atoms,
        std::vector<TopologyResidue>& residues,
        std::vector<TopologySegment>& segments,
        std::unordered_map<std::string, int>& segment_map,
        std::map<std::tuple<std::string, int, std::string>, int>& residue_map,
        std::map<std::tuple<std::string, int, std::string, std::string>, int>& atom_map,
        const std::string& name, const std::string& type, double charge, double mass,
        const std::string& residue_name, int residue_number, const std::string& segment_name
    ) {
        // Ensure we have the segment
        int segment_id;
        auto segment_it = segment_map.find(segment_name);
        if (segment_it == segment_map.end()) {
            TopologySegment segment;
            segment.id = static_cast<int>(segments.size());
            segment.name = segment_name;
            segment_id = segment.id;
            segments.push_back(segment);
            segment_map[segment_name] = segment_id;
        } else {
            segment_id = segment_it->second;
        }

        // Ensure we have the residue
        int residue_id;
        auto residue_key = std::make_tuple(residue_name, residue_number, segment_name);
        auto residue_it = residue_map.find(residue_key);
        if (residue_it == residue_map.end()) {
            TopologyResidue residue;
            residue.id = static_cast<int>(residues.size());
            residue.name = residue_name;
            residue.number = residue_number;
            residue.segment = segment_name;
            residue_id = residue.id;
            residues.push_back(residue);
            residue_map[residue_key] = residue_id;
            segments[segment_id].residues.push_back(residue_id);
        } else {
            residue_id = residue_it->second;
        }

        // Create and add the atom
        TopologyAtom atom;
        atom.id = static_cast<int>(atoms.size());
        atom.name = name;
        atom.type = type;
        atom.charge = charge;
        atom.mass = mass;
        atom.residue_id = residue_id;
        atom.segment_id = segment_id;

        int atom_id = static_cast<int>(atoms.size());
        atoms.push_back(atom);
        atom_map[std::make_tuple(residue_name, residue_number, segment_name, name)] = atom_id;
        residues[residue_id].atoms.push_back(atom_id);

        return atom_id;
    }

    // Add connectivity operations
    static inline void addBond(std::vector<TopologyBond>& bonds, 
                              const std::vector<TopologyAtom>& atoms,
                              int atom1, int atom2, double length = 0.0, 
                              double force_constant = 0.0, int function_type = 1) {
        if (!hasAtom(atoms, atom1) || !hasAtom(atoms, atom2)) {
            throw std::out_of_range("Invalid atom indices in add_bond");
        }
        TopologyBond bond;
        bond.atom1 = atom1;
        bond.atom2 = atom2;
        bond.length = length;
        bond.force_constant = force_constant;
        bond.function_type = function_type;
        bonds.push_back(bond);
    }

    static inline void addAngle(std::vector<TopologyAngle>& angles,
                               const std::vector<TopologyAtom>& atoms,
                               int atom1, int atom2, int atom3, double angle = 0.0,
                               double force_constant = 0.0, int function_type = 1) {
        if (!hasAtom(atoms, atom1) || !hasAtom(atoms, atom2) || !hasAtom(atoms, atom3)) {
            throw std::out_of_range("Invalid atom indices in add_angle");
        }
        TopologyAngle ang;
        ang.atom1 = atom1;
        ang.atom2 = atom2;
        ang.atom3 = atom3;
        ang.angle = angle;
        ang.force_constant = force_constant;
        ang.function_type = function_type;
        angles.push_back(ang);
    }

    static inline void addDihedral(std::vector<TopologyDihedral>& dihedrals,
                                  const std::vector<TopologyAtom>& atoms,
                                  int atom1, int atom2, int atom3, int atom4, int multiplicity = 1,
                                  double angle = 0.0, double force_constant = 0.0, 
                                  bool improper = false, int function_type = 1) {
        if (!hasAtom(atoms, atom1) || !hasAtom(atoms, atom2) || !hasAtom(atoms, atom3) || !hasAtom(atoms, atom4)) {
            throw std::out_of_range("Invalid atom indices in add_dihedral");
        }
        TopologyDihedral dihedral;
        dihedral.atom1 = atom1;
        dihedral.atom2 = atom2;
        dihedral.atom3 = atom3;
        dihedral.atom4 = atom4;
        dihedral.multiplicity = multiplicity;
        dihedral.angle = angle;
        dihedral.force_constant = force_constant;
        dihedral.improper = improper;
        dihedral.function_type = function_type;
        dihedrals.push_back(dihedral);
    }

    // Add special features
    static inline void addDonor(std::vector<TopologyDonor>& donors,
                               const std::vector<TopologyAtom>& atoms,
                               int donor, int hydrogen) {
        if (!hasAtom(atoms, donor) || !hasAtom(atoms, hydrogen)) {
            throw std::out_of_range("Invalid atom indices in add_donor");
        }
        TopologyDonor d;
        d.donor_atom = donor;
        d.hydrogen_atom = hydrogen;
        donors.push_back(d);
    }

    static inline void addAcceptor(std::vector<TopologyAcceptor>& acceptors,
                                  const std::vector<TopologyAtom>& atoms,
                                  int acceptor) {
        if (!hasAtom(atoms, acceptor)) {
            throw std::out_of_range("Invalid atom index in add_acceptor");
        }
        TopologyAcceptor a;
        a.acceptor_atom = acceptor;
        acceptors.push_back(a);
    }

    static inline void addGroup(std::vector<TopologyGroup>& groups,
                               const std::vector<TopologyAtom>& atoms,
                               int id, const std::vector<int>& atoms_list, const std::string& type = "") {
        for (int atom : atoms_list) {
            if (!hasAtom(atoms, atom)) {
                throw std::out_of_range("Invalid atom index in add_group");
            }
        }
        TopologyGroup group;
        group.id = id;
        group.atoms = atoms_list;
        group.type = type;
        groups.push_back(group);
    }

    static inline void addCmap(std::vector<TopologyCmap>& cmaps,
                              const std::vector<TopologyAtom>& atoms,
                              const std::array<int, 8>& atoms_array) {
        for (int atom : atoms_array) {
            if (atom >= 0 && !hasAtom(atoms, atom)) {
                throw std::out_of_range("Invalid atom index in add_cmap");
            }
        }
        TopologyCmap cmap;
        cmap.atoms = atoms_array;
        cmap.function_type = 1;
        cmaps.push_back(cmap);
    }

    static inline void addCmap(std::vector<TopologyCmap>& cmaps,
                              const std::vector<TopologyAtom>& /*atoms*/,
                              const std::array<int, 5>& atoms_array, int function_type = 1) {
        std::array<int, 8> charmm_atoms;
        for (int i = 0; i < 5; ++i) {
            charmm_atoms[i] = atoms_array[i];
        }
        for (int i = 5; i < 8; ++i) {
            charmm_atoms[i] = -1;
        }
        TopologyCmap cmap;
        cmap.atoms = charmm_atoms;
        cmap.function_type = function_type;
        cmaps.push_back(cmap);
    }


    static inline void addExclusion(std::map<int, std::set<int>>& exclusions,
                                   const std::vector<TopologyAtom>& atoms,
                                   int atom1, int atom2) {
        if (!hasAtom(atoms, atom1) || !hasAtom(atoms, atom2)) {
            throw std::out_of_range("Invalid atom indices in add_nonbonded_exclusion");
        }
        exclusions[atom1].insert(atom2);
        exclusions[atom2].insert(atom1);
    }

    // Helper operations
    static inline int addSegment(std::vector<TopologySegment>& segments,
                                std::unordered_map<std::string, int>& segment_map,
                                const std::string& name) {
        TopologySegment segment;
        segment.id = segments.size();
        segment.name = name;

        int segment_id = segments.size();
        segments.push_back(segment);
        segment_map[name] = segment_id;

        return segment_id;
    }

    static inline int addResidue(std::vector<TopologyResidue>& residues,
                                std::vector<TopologySegment>& segments,
                                std::unordered_map<std::string, int>& segment_map,
                                std::map<std::tuple<std::string, int, std::string>, int>& residue_map,
                                const std::string& name, int number, const std::string& segment) {
        // Ensure we have the segment
        int segment_id;
        auto segment_it = segment_map.find(segment);
        if (segment_it == segment_map.end()) {
            segment_id = addSegment(segments, segment_map, segment);
        } else {
            segment_id = segment_it->second;
        }

        // Create and add the residue
        TopologyResidue residue;
        residue.id = residues.size();
        residue.name = name;
        residue.number = number;
        residue.segment = segment;

        int residue_id = residues.size();
        residues.push_back(residue);
        residue_map[std::make_tuple(name, number, segment)] = residue_id;
        segments[segment_id].residues.push_back(residue_id);

        return residue_id;
    }

private:
    static inline bool hasAtom(const std::vector<TopologyAtom>& atoms, int index) {
        return index >= 0 && index < static_cast<int>(atoms.size());
    }
};

} // namespace topology
} // namespace model
} // namespace pygcmc