#pragma once

#ifndef PYGCMC_MODEL_ATOM_MAIN_HPP
#define PYGCMC_MODEL_ATOM_MAIN_HPP

#include "AtomCore.hpp"
#include <string>
#include <sstream>
#include <iomanip>

namespace pygcmc {
namespace model {
namespace atom {

/**
 * @brief Extended Atom class with backward compatibility and utility functions
 * This class extends AtomCore with additional functionality and maintains compatibility
 */
class Atom : public AtomCore {
public:
    // Inherit constructors
    using AtomCore::AtomCore;

    // Additional getters for backward compatibility
    const std::string& get_element() const noexcept { return element; }
    const std::string& get_charge_string() const noexcept { return chargestr; }

    // Additional setters for PDB compatibility
    void set_occupancy(double occ) {
        if (occ < 0.0 || occ > 1.0) {
            throw std::invalid_argument("Occupancy must be between 0 and 1");
        }
        occupancy = occ;
    }

    void set_tempfactor(double temp) {
        if (!std::isfinite(temp)) {
            throw std::invalid_argument("Invalid temperature factor");
        }
        tempfactor = temp;
    }

    void set_element(const std::string& elem) {
        element = elem;
    }

    void set_charge_string(const std::string& chg) {
        chargestr = chg;
    }

    void set_altloc(char alt) { altloc = alt; }
    void set_inscode(char ins) { inscode = ins; }

    // PDB format utilities
    static std::string format_pdb_atom_name(const std::string& name) {
        // Left align atom name according to PDB format
        // Element symbols are right-justified in columns 13-14
        if (name.length() >= 4) return name;
        
        // Check if first character is a digit (indicating branch)
        if (!name.empty() && std::isdigit(name[0])) {
            return name;  // Left justify if starts with digit
        }
        
        // Right justify element symbol
        std::string result(4, ' ');
        if (name.length() == 1) {
            result[1] = name[0];  // Single character element
        } else if (name.length() > 1) {
            result[0] = name[0];  // Two character element
            result[1] = name[1];
        }
        
        // Add remaining characters
        for (size_t i = 2; i < name.length() && i < 4; ++i) {
            result[i] = name[i];
        }
        
        return result;
    }

    std::string get_formatted_atom_name() const {
        return format_pdb_atom_name(type);
    }

    std::string get_residue_id() const {
        // Combine residue number and insertion code (e.g., "153A")
        if (inscode == ' ') {
            return std::to_string(ires);
        }
        return std::to_string(ires) + inscode;
    }

    void set_residue_id(const std::string& resid) {
        // Parse residue ID (e.g., "153A" -> ires=153, inscode='A')
        size_t numLen = 0;
        try {
            ires = std::stoi(resid, &numLen);
        } catch (const std::exception&) {
            throw std::invalid_argument("Invalid residue ID format");
        }
        
        if (numLen < resid.length()) {
            inscode = resid[numLen];
        } else {
            inscode = ' ';
        }
    }

    // Advanced utility methods
    std::string get_atom_identifier() const {
        // Generate unique atom identifier: segid:resname:ires:type
        std::stringstream ss;
        ss << segid << ":" << resname << ":" << ires << ":" << type;
        return ss.str();
    }

    std::string get_pdb_record() const {
        // Generate PDB ATOM/HETATM record
        std::stringstream ss;
        
        // Record type
        ss << (hetatm ? "HETATM" : "ATOM  ");
        
        // Atom serial number (5 chars, right-aligned)
        ss << std::setw(5) << std::right << bynu;
        
        // Space
        ss << " ";
        
        // Atom name (4 chars, formatted)
        ss << get_formatted_atom_name();
        
        // Alternate location (1 char)
        ss << altloc;
        
        // Residue name (3 chars, left-aligned)
        ss << std::setw(3) << std::left << resname;
        
        // Space + Chain ID
        ss << " " << chain;
        
        // Residue sequence number (4 chars, right-aligned)
        ss << std::setw(4) << std::right << ires;
        
        // Insertion code
        ss << inscode;
        
        // Spaces (3 chars)
        ss << "   ";
        
        // Coordinates (8.3 format each)
        ss << std::fixed << std::setprecision(3);
        ss << std::setw(8) << std::right << coor[0];
        ss << std::setw(8) << std::right << coor[1];
        ss << std::setw(8) << std::right << coor[2];
        
        // Occupancy and temperature factor
        ss << std::setw(6) << std::setprecision(2) << occupancy;
        ss << std::setw(6) << std::setprecision(2) << tempfactor;
        
        // Spaces (10 chars)
        ss << "          ";
        
        // Element symbol (2 chars, right-aligned)
        if (!element.empty()) {
            ss << std::setw(2) << std::right << element;
        } else {
            ss << "  ";
        }
        
        // Charge (2 chars)
        if (!chargestr.empty()) {
            ss << std::setw(2) << std::right << chargestr;
        }
        
        return ss.str();
    }

    // Distance calculation utilities
    double distance_to(const Atom& other) const {
        double dx = coor[0] - other.coor[0];
        double dy = coor[1] - other.coor[1];
        double dz = coor[2] - other.coor[2];
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }

    double distance_squared_to(const Atom& other) const {
        double dx = coor[0] - other.coor[0];
        double dy = coor[1] - other.coor[1];
        double dz = coor[2] - other.coor[2];
        return dx*dx + dy*dy + dz*dz;
    }

    // Comparison operators for sorting/searching
    bool operator<(const Atom& other) const {
        if (segid != other.segid) return segid < other.segid;
        if (ires != other.ires) return ires < other.ires;
        return bynu < other.bynu;
    }

    bool operator==(const Atom& other) const {
        return bynu == other.bynu && 
               type == other.type && 
               resname == other.resname &&
               ires == other.ires && 
               segid == other.segid &&
               std::abs(coor[0] - other.coor[0]) < 1e-9 &&
               std::abs(coor[1] - other.coor[1]) < 1e-9 &&
               std::abs(coor[2] - other.coor[2]) < 1e-9;
    }

    bool operator!=(const Atom& other) const {
        return !(*this == other);
    }


};

// Utility functions for atom collections
namespace utils {

/**
 * @brief Find atoms by type in a collection
 */
template<typename Container>
auto find_atoms_by_type(const Container& atoms, const std::string& type) {
    std::vector<typename Container::value_type> result;
    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(result),
                [&type](const auto& atom) { return atom.get_type() == type; });
    return result;
}

/**
 * @brief Find atoms by residue in a collection
 */
template<typename Container>
auto find_atoms_by_residue(const Container& atoms, const std::string& resname, int ires) {
    std::vector<typename Container::value_type> result;
    std::copy_if(atoms.begin(), atoms.end(), std::back_inserter(result),
                [&resname, ires](const auto& atom) { 
                    return atom.get_resname() == resname && atom.get_ires() == ires; 
                });
    return result;
}

} // namespace utils

/**
 * @brief Validate an atom
 * @param atom The atom to validate
 * @return true if the atom is valid, false otherwise
 */
inline bool validate_atom(const Atom& atom) {
    return atom.is_valid();
}

} // namespace atom
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_ATOM_MAIN_HPP 