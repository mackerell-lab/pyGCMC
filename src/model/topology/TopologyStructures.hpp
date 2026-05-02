#pragma once

#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <array>
#include <optional>
#include <tuple>

namespace pygcmc {
namespace model {
namespace topology {

/**
 * @brief Represents an atom in the topology
 */
struct TopologyAtom {
    int id;                 ///< Atom ID (0-based)
    std::string name;       ///< Atom name (e.g., "CA", "N", "O")
    std::string type;       ///< Atom type (e.g., "CT1", "NH1")
    double charge;          ///< Atomic charge
    double mass;            ///< Atomic mass
    int residue_id;        ///< ID of the residue this atom belongs to
    int segment_id;        ///< ID of the segment this atom belongs to
    double alpha = 0.0;     ///< Polarizability (Drude force field)
    double thole = 0.0;     ///< Thole screening parameter (Drude force field)

    // Getter methods for consistency with AtomCore interface
    double get_alpha() const { return alpha; }
    double get_thole() const { return thole; }

    // Setter methods for Drude parameters
    void set_alpha(double a) { alpha = a; }
    void set_thole(double t) { thole = t; }
    void set_drude_params(double a, double t) { alpha = a; thole = t; }
};

/**
 * @brief Represents a residue in the topology
 */
struct TopologyResidue {
    int id;                 ///< Residue ID (0-based)
    std::string name;       ///< Residue name (e.g., "ALA", "GLY")
    int number;             ///< Residue number from PDB/PSF
    std::vector<int> atoms; ///< Indices of atoms in this residue
    std::string segment;    ///< Segment name this residue belongs to
};

/**
 * @brief Represents a segment in the topology
 */
struct TopologySegment {
    int id;                 ///< Segment ID (0-based)
    std::string name;       ///< Segment name (e.g., "PROT", "MEMB")
    std::vector<int> residues; ///< Indices of residues in this segment
};

/**
 * @brief Represents a bond between two atoms
 */
struct TopologyBond {
    int atom1;              ///< Index of first atom
    int atom2;              ///< Index of second atom
    double length;          ///< Equilibrium bond length (optional)
    double force_constant;  ///< Bond force constant (optional)
    int function_type = 1;    // Default to GROMACS function type 1
};

/**
 * @brief Represents an angle between three atoms
 */
struct TopologyAngle {
    int atom1;              ///< Index of first atom
    int atom2;              ///< Index of central atom
    int atom3;              ///< Index of third atom
    double angle;           ///< Equilibrium angle in degrees (optional)
    double force_constant;  ///< Angle force constant (optional)
    int function_type = 1;   // Default to GROMACS function type 1
};

/**
 * @brief Represents a dihedral angle between four atoms
 */
struct TopologyDihedral {
    int atom1;              ///< Index of first atom
    int atom2;              ///< Index of second atom
    int atom3;              ///< Index of third atom
    int atom4;              ///< Index of fourth atom
    int multiplicity;       ///< Dihedral multiplicity
    double angle;           ///< Equilibrium angle in degrees
    double force_constant;  ///< Dihedral force constant
    bool improper;          ///< Whether this is an improper dihedral
    int function_type = 1;   // Default to GROMACS function type 1
};

/**
 * @brief Represents a hydrogen bond donor
 */
struct TopologyDonor {
    int donor_atom;     ///< Index of the donor atom
    int hydrogen_atom;  ///< Index of the hydrogen atom
};

/**
 * @brief Represents a hydrogen bond acceptor
 */
struct TopologyAcceptor {
    int acceptor_atom;  ///< Index of the acceptor atom
};

/**
 * @brief Represents a group in the topology
 */
struct TopologyGroup {
    int id;                    ///< Group ID
    std::vector<int> atoms;    ///< Indices of atoms in this group
    std::string type;          ///< Group type (e.g., "WATER")
};

/**
 * @brief Represents a CMAP (correction map) term
 */
struct TopologyCmap {
    std::array<int, 8> atoms;  ///< 8 atoms involved in CHARMM CMAP term
    int function_type = 1;     ///< CMAP function type (default: 1)
};

} // namespace topology
} // namespace model
} // namespace pygcmc
