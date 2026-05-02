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

    // Convenience method for get_name() to match common API expectations
    std::string get_name() const {
        return type;
    }

    // Additional convenience properties for enhanced test compatibility
    std::string atom_name() const {
        return type;
    }

    std::string residue_name() const {
        return resname;
    }

    int residue_number() const {
        return ires;
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


} // namespace atom
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_ATOM_MAIN_HPP
