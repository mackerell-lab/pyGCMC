#pragma once

#ifndef PYGCMC_MODEL_ATOM_MAIN_HPP
#define PYGCMC_MODEL_ATOM_MAIN_HPP

#include "AtomCore.hpp"
#include <sstream>
#include <iomanip>

namespace pygcmc {
namespace model {

/**
 * @brief Extended Atom class with PDB formatting and utility functions
 * @details Extends AtomCore with additional functionality while maintaining compatibility
 */
class Atom : public AtomCore, public ICloneable<Atom>, public ISerializable {
public:
    // Inherit constructors
    using AtomCore::AtomCore;

    // Default constructor
    Atom() = default;

    // Copy constructor and assignment
    Atom(const Atom&) = default;
    Atom& operator=(const Atom&) = default;
    Atom(Atom&&) = default;
    Atom& operator=(Atom&&) = default;

    // Constructor from AtomCore
    explicit Atom(const AtomCore& core) : AtomCore(core) {}

    // Additional setters for PDB compatibility
    void set_occupancy(double occ) {
        if (occ < 0.0 || occ > 1.0) {
            throw std::invalid_argument("Occupancy must be between 0 and 1");
        }
        occupancy_ = occ;
    }

    void set_tempfactor(double temp) {
        if (!std::isfinite(temp)) {
            throw std::invalid_argument("Invalid temperature factor");
        }
        tempfactor_ = temp;
    }

    void set_charge_string(const std::string& chg) {
        chargestr_ = chg;
    }

    // PDB format utilities
    static std::string format_pdb_atom_name(const std::string& name) {
        return utils::format::format_pdb_atom_name(name);
    }

    std::string get_formatted_atom_name() const {
        return format_pdb_atom_name(type_);
    }

    std::string get_residue_id() const {
        if (inscode_ == ' ') {
            return std::to_string(ires_);
        }
        return std::to_string(ires_) + inscode_;
    }

    void set_residue_id(const std::string& resid) {
        size_t numLen = 0;
        try {
            ires_ = std::stoi(resid, &numLen);
        } catch (const std::exception&) {
            throw std::invalid_argument("Invalid residue ID format");
        }
        
        if (numLen < resid.length()) {
            inscode_ = resid[numLen];
        } else {
            inscode_ = ' ';
        }
    }

    // Distance calculations
    double distance_to(const Atom& other) const {
        return utils::math::distance(coor_, other.get_coor());
    }

    double distance_squared_to(const Atom& other) const {
        return utils::math::distance_squared(coor_, other.get_coor());
    }

    // Comparison operators
    bool operator==(const Atom& other) const {
        return bynu_ == other.bynu_ && 
               type_ == other.type_ && 
               resname_ == other.resname_ &&
               ires_ == other.ires_ && 
               segid_ == other.segid_ &&
               utils::compare::coordinates_equal(coor_, other.coor_);
    }

    bool operator!=(const Atom& other) const {
        return !(*this == other);
    }

    // Hash support
    std::size_t hash() const {
        return utils::hash::combine_hash(
            bynu_, utils::hash::string_hash(type_), 
            utils::hash::string_hash(resname_), ires_,
            utils::hash::string_hash(segid_),
            utils::hash::coordinate_hash(coor_)
        );
    }

    // ICloneable interface
    std::unique_ptr<Atom> clone() const override {
        return std::make_unique<Atom>(*this);
    }

    // ISerializable interface
    std::string serialize() const override {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3);
        oss << "ATOM:" << bynu_ << ":" << type_ << ":" << resname_ 
            << ":" << ires_ << ":" << segid_ << ":" << chain_
            << ":" << coor_[0] << ":" << coor_[1] << ":" << coor_[2]
            << ":" << mass_ << ":" << charge_;
        return oss.str();
    }

    bool deserialize(const std::string& data) override {
        std::istringstream iss(data);
        std::string token;
        
        if (!std::getline(iss, token, ':') || token != "ATOM") return false;
        if (!std::getline(iss, token, ':')) return false;
        bynu_ = std::stoi(token);
        if (!std::getline(iss, token, ':')) return false;
        type_ = token;
        if (!std::getline(iss, resname_, ':')) return false;
        if (!std::getline(iss, token, ':')) return false;
        ires_ = std::stoi(token);
        if (!std::getline(iss, segid_, ':')) return false;
        if (!std::getline(iss, token, ':')) return false;
        chain_ = token.empty() ? ' ' : token[0];
        if (!std::getline(iss, token, ':')) return false;
        coor_[0] = std::stod(token);
        if (!std::getline(iss, token, ':')) return false;
        coor_[1] = std::stod(token);
        if (!std::getline(iss, token, ':')) return false;
        coor_[2] = std::stod(token);
        if (!std::getline(iss, token, ':')) return false;
        mass_ = std::stod(token);
        if (!std::getline(iss, token)) return false;
        charge_ = std::stod(token);
        
        return true;
    }

    std::string get_type_name() const override {
        return "Atom";
    }

    // PDB format output
    std::string to_pdb_string() const {
        std::ostringstream oss;
        oss << std::left;
        oss << std::setw(6) << (hetatm_ ? "HETATM" : "ATOM");
        oss << std::right << std::setw(5) << bynu_;
        oss << " ";
        oss << std::left << std::setw(4) << get_formatted_atom_name();
        oss << std::setw(1) << altloc_;
        oss << std::setw(3) << resname_;
        oss << " ";
        oss << std::setw(1) << chain_;
        oss << std::right << std::setw(4) << ires_;
        oss << std::setw(1) << inscode_;
        oss << "   ";
        oss << std::fixed << std::setprecision(3);
        oss << std::right << std::setw(8) << coor_[0];
        oss << std::setw(8) << coor_[1];
        oss << std::setw(8) << coor_[2];
        oss << std::setprecision(2);
        oss << std::setw(6) << occupancy_;
        oss << std::setw(6) << tempfactor_;
        oss << "          ";
        oss << std::left << std::setw(2) << element_;
        oss << std::setw(2) << chargestr_;
        
        return oss.str();
    }

    // Enhanced validation with detailed error messages
    std::string get_validation_error() const override {
        if (bynu_ <= 0) return "Invalid atom number";
        if (!utils::validate::is_valid_name(type_)) return "Invalid atom type";
        if (!utils::validate::is_valid_name(resname_)) return "Invalid residue name";
        if (ires_ <= 0) return "Invalid residue number";
        if (!utils::validate::is_valid_mass(mass_)) return "Invalid mass";
        if (!utils::validate::is_valid_charge(charge_)) return "Invalid charge";
        if (!std::all_of(coor_.begin(), coor_.end(), utils::validate::is_valid_coordinate)) {
            return "Invalid coordinates";
        }
        if (occupancy_ < 0.0 || occupancy_ > 1.0) return "Invalid occupancy";
        if (!std::isfinite(tempfactor_)) return "Invalid temperature factor";
        
        return "";
    }
};

} // namespace model
} // namespace pygcmc

// Hash specialization for std::unordered_map support
namespace std {
    template<>
    struct hash<pygcmc::model::Atom> {
        std::size_t operator()(const pygcmc::model::Atom& atom) const {
            return atom.hash();
        }
    };
}

#endif // PYGCMC_MODEL_ATOM_MAIN_HPP 