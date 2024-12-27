// modules/core/include/pygcmc/core/io/parser_common.hpp

#ifndef PYGCMC_CORE_IO_PARSER_COMMON_HPP
#define PYGCMC_CORE_IO_PARSER_COMMON_HPP

#include <vector>
#include <array>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <utility>
#include <stdexcept>
#include <functional>
#include <cmath>
#include <tuple>
#include "pygcmc/core/utils.hpp"

namespace pygcmc {
namespace core {
namespace io {

/**
 * @brief Base parser exception class
 */
class ParserError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/**
 * @brief File format error
 */
class FormatError : public ParserError {
public:
    using ParserError::ParserError;
};

/**
 * @brief File access error
 */
class FileError : public ParserError {
public:
    using ParserError::ParserError;
};

/**
 * @brief Common atom structure for PDB files
 */
struct PDBAtom {
    int serial;              ///< Atom serial number
    std::string name;        ///< Atom name
    std::string residue;     ///< Residue name
    int sequence;            ///< Residue sequence number
    char chain;              ///< Chain identifier
    char alt_loc;            ///< Alternate location indicator
    char insertion_code;     ///< Insertion code
    double x, y, z;          ///< Atomic coordinates
    double occupancy;        ///< Occupancy
    double temp_factor;      ///< Temperature factor
    std::string element;     ///< Element symbol
    std::string charge;      ///< Charge
    std::string type;        ///< Atom type

    PDBAtom(int serial_ = 0, const std::string& name_ = "", 
            const std::string& residue_ = "", int sequence_ = 0,
            char chain_ = ' ', char alt_loc_ = ' ', char insertion_code_ = ' ',
            double x_ = 0.0, double y_ = 0.0, double z_ = 0.0,
            double occupancy_ = 0.0, double temp_factor_ = 0.0,
            const std::string& element_ = "", const std::string& charge_ = "",
            const std::string& type_ = "")
        : serial(serial_), name(name_), residue(residue_),
          sequence(sequence_), chain(chain_), alt_loc(alt_loc_), insertion_code(insertion_code_),
          x(x_), y(y_), z(z_), occupancy(occupancy_), temp_factor(temp_factor_),
          element(element_), charge(charge_), type(type_) {}
    
    bool is_valid() const {
        return serial > 0 && !name.empty() && !residue.empty() &&
               std::isfinite(x) && std::isfinite(y) && std::isfinite(z) &&
               std::isfinite(occupancy) && std::isfinite(temp_factor) &&
               !element.empty() && !type.empty();
    }

    std::array<double, 3> position() const {
        return {x, y, z};
    }
};


/**
 * @brief ITP Atom structure
 */
struct ITPAtom {
    std::string name;    ///< Atom name
    std::string type;    ///< Atom type
    int resid;           ///< Residue ID
    std::string resname; ///< Residue name
    double charge;       ///< Atomic charge

    ITPAtom(const std::string& name_ = "", const std::string& type_ = "",
            int resid_ = 0, const std::string& resname_ = "", double charge_ = 0.0)
        : name(name_), type(type_), resid(resid_), resname(resname_), charge(charge_) {}

    bool is_valid() const {
        return !name.empty() && !type.empty() && resid > 0 && !resname.empty() && std::isfinite(charge);
    }
};

/**
 * @brief PSF Bond structure
 */
struct PSFBond {
    int atom1; ///< First atom index
    int atom2; ///< Second atom index

    PSFBond(int a1 = 0, int a2 = 0) : atom1(a1), atom2(a2) {}

    bool is_valid() const {
        return atom1 > 0 && atom2 > 0 && atom1 != atom2;
    }
};

/**
 * @brief PSF Atom structure
 */
struct PSFAtom {
    int id;              ///< Atom ID
    int residue_id;      ///< Residue ID
    std::string name;    ///< Atom name
    std::string type;    ///< Atom type
    double charge;       ///< Atom charge
    double mass;         ///< Atom mass

    PSFAtom(int id_ = 0, int residue_id_ = 0, const std::string& name_ = "", 
            const std::string& type_ = "", double charge_ = 0.0, double mass_ = 0.0)
        : id(id_), residue_id(residue_id_), name(name_), type(type_), charge(charge_), mass(mass_) {}

    bool is_valid() const {
        return id > 0 && !name.empty() && !type.empty() &&
               std::isfinite(charge) && std::isfinite(mass);
    }
};

/**
 * @brief PSF Topology structure
 */
struct PSFTopology {
    std::vector<PSFAtom> atoms;  ///< List of atoms in the PSF
    std::vector<PSFBond> bonds;  ///< List of bonds in the PSF

    bool is_valid() const {
        // Check if we have any atoms and bonds
        if (atoms.empty()) {
            return false;
        }

        // Validate all atoms
        for (const auto& atom : atoms) {
            if (!atom.is_valid()) {
                return false;
            }
        }

        // Create a set of valid atom IDs
        std::unordered_set<int> valid_atom_ids;
        for (const auto& atom : atoms) {
            valid_atom_ids.insert(atom.id);
        }

        // Check bonds
        std::set<std::pair<int, int>> bond_pairs;
        for (const auto& bond : bonds) {
            // Check if bond is valid
            if (!bond.is_valid()) {
                return false;
            }

            // Check if bond references valid atoms
            if (valid_atom_ids.find(bond.atom1) == valid_atom_ids.end() ||
                valid_atom_ids.find(bond.atom2) == valid_atom_ids.end()) {
                return false;
            }

            // Check for duplicate bonds
            int min_atom = std::min(bond.atom1, bond.atom2);
            int max_atom = std::max(bond.atom1, bond.atom2);
            auto bond_pair = std::make_pair(min_atom, max_atom);
            
            if (!bond_pairs.insert(bond_pair).second) {
                return false;  // Duplicate bond found
            }
        }

        return true;
    }
};

/**
 * @brief Topology Atom Type structure
 */
struct TopAtomType {
    std::string name; ///< Atom type name
    int type;         ///< Atom type number
    double charge;    ///< Atomic charge
    double mass;      ///< Atomic mass

    TopAtomType(const std::string& name_ = "", int type_ = 0,
               double charge_ = 0.0, double mass_ = 0.0)
        : name(name_), type(type_), charge(charge_), mass(mass_) {}

    bool is_valid() const {
        return !name.empty() && type > 0 && std::isfinite(charge) && std::isfinite(mass);
    }
};

/**
 * @brief Force field parameter pair structure for non-bonded interactions
 */
struct ForceFieldPair {
    double param1; ///< First parameter (typically sigma in nm)
    double param2; ///< Second parameter (typically epsilon in kJ/mol)

    ForceFieldPair(double p1 = 0.0, double p2 = 0.0) 
        : param1(p1), param2(p2) {}

    bool is_valid() const {
        return std::isfinite(param1) && std::isfinite(param2) && 
               param1 >= 0.0;  // sigma should be non-negative
    }

    /**
     * @brief Combines two ForceFieldPairs using arithmetic mixing rules
     * @param other The other ForceFieldPair to combine with
     * @return A new ForceFieldPair with combined parameters
     */
    ForceFieldPair combine_arithmetic(const ForceFieldPair& other) const {
        return ForceFieldPair(
            0.5 * (param1 + other.param1),     // arithmetic mean of sigma
            std::sqrt(param2 * other.param2)   // geometric mean of epsilon
        );
    }

    /**
     * @brief Combines two ForceFieldPairs using geometric mixing rules
     * @param other The other ForceFieldPair to combine with
     * @return A new ForceFieldPair with combined parameters
     */
    ForceFieldPair combine_geometric(const ForceFieldPair& other) const {
        return ForceFieldPair(
            std::sqrt(param1 * other.param1),   // geometric mean of sigma
            std::sqrt(param2 * other.param2)    // geometric mean of epsilon
        );
    }

    /**
     * @brief Scales the parameters by a factor
     * @param factor The scaling factor
     * @return A new ForceFieldPair with scaled parameters
     */
    ForceFieldPair scale(double factor) const {
        return ForceFieldPair(param1 * factor, param2 * factor);
    }
};

struct PairStringHash {
    std::size_t operator()(const std::pair<std::string, std::string>& p) const {
        return std::hash<std::string>()(p.first) ^ (std::hash<std::string>()(p.second) << 1);
    }
};

using NBMap = std::unordered_map<std::string, ForceFieldPair>;
using NBFixMap = std::unordered_map<std::pair<std::string, std::string>, ForceFieldPair, PairStringHash>;

struct Topology {
    std::vector<TopAtomType> atom_types;

    bool is_valid() const {
        if (atom_types.empty()) return false;
        for(const auto& atom_type : atom_types){
            if(!atom_type.is_valid()) return false;
        }
        return true;
    }
};

} // namespace io
} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_IO_PARSER_COMMON_HPP
