#pragma once

#ifndef PYGCMC_MODEL_ATOM_CORE_HPP
#define PYGCMC_MODEL_ATOM_CORE_HPP

#include <string>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include "../common/ModelInterface.hpp"
#include "../common/ModelConstants.hpp"
#include "../common/ModelUtils.hpp"

namespace pygcmc {
namespace model {

/**
 * @brief Core Atom class with essential data and operations
 * @details Lightweight atom representation with CHARMM naming conventions
 */
class AtomCore : public IValidatable, public IIdentifiable {
public:
    // Default constructor
    AtomCore() : 
        bynu_(0), ires_(0), iseg_(0), igro_(0),
        altloc_(' '), chain_(' '), inscode_(' '),
        coor_{0.0, 0.0, 0.0},
        occupancy_(1.0), tempfactor_(0.0),
        wmain_(1.0), wcomp_(1.0), mass_(0.0), charge_(0.0),
        radius_(0.0), alpha_(0.0),
        eps_(constants::INVALID_VALUE), rmin_(constants::INVALID_VALUE),
        fbeta_(0.0), move_(1), ignore_(0), constrain_(0),
        hetatm_(false), initial_(false) {}

    // Parameterized constructor
    AtomCore(int bynu, const std::string& type, const std::string& resname,
             int ires, const std::string& segid = "", int iseg = 0,
             double x = 0.0, double y = 0.0, double z = 0.0,
             double wmain = 1.0, double mass = 0.0, double charge = 0.0,
             const std::string& chem = "", bool hetatm = false) :
        bynu_(bynu), type_(type), resname_(resname), ires_(ires),
        segid_(segid), iseg_(iseg), igro_(0),
        altloc_(' '), chain_(' '), inscode_(' '),
        coor_{x, y, z}, wmain_(wmain), mass_(mass), charge_(charge),
        chem_(chem), eps_(constants::INVALID_VALUE), rmin_(constants::INVALID_VALUE),
        hetatm_(hetatm) {}

    // CHARMM standard getters
    int get_bynu() const noexcept { return bynu_; }
    const std::string& get_type() const noexcept { return type_; }
    const std::string& get_resname() const noexcept { return resname_; }
    int get_ires() const noexcept { return ires_; }
    const std::string& get_segid() const noexcept { return segid_; }
    int get_iseg() const noexcept { return iseg_; }
    int get_igro() const noexcept { return igro_; }
    const std::string& get_chem() const noexcept { return chem_; }
    double get_wmain() const noexcept { return wmain_; }
    double get_mass() const noexcept { return mass_; }
    double get_charge() const noexcept { return charge_; }
    
    // Coordinate access
    const std::array<double, 3>& get_coor() const noexcept { return coor_; }
    double get_x() const noexcept { return coor_[0]; }
    double get_y() const noexcept { return coor_[1]; }
    double get_z() const noexcept { return coor_[2]; }

    // Force field parameters
    double get_eps() const noexcept { return eps_; }
    double get_rmin() const noexcept { return rmin_; }

    // PDB specific getters
    char get_altloc() const noexcept { return altloc_; }
    char get_inscode() const noexcept { return inscode_; }
    char get_chain() const noexcept { return chain_; }
    bool is_hetatm() const noexcept { return hetatm_; }
    double get_occupancy() const noexcept { return occupancy_; }
    double get_tempfactor() const noexcept { return tempfactor_; }
    const std::string& get_element() const noexcept { return element_; }
    const std::string& get_charge_string() const noexcept { return chargestr_; }

    // Core setters with validation
    void set_coor(double x, double y, double z) {
        if (!utils::validate::is_valid_coordinate(x) || 
            !utils::validate::is_valid_coordinate(y) || 
            !utils::validate::is_valid_coordinate(z)) {
            throw std::invalid_argument("Invalid coordinates");
        }
        coor_ = {x, y, z};
    }

    void set_mass_charge(double m, double q) {
        if (!utils::validate::is_valid_mass(m) || !utils::validate::is_valid_charge(q)) {
            throw std::invalid_argument("Invalid mass or charge");
        }
        mass_ = m;
        charge_ = q;
    }

    void set_lj_params(double epsilon, double r) {
        if (!std::isfinite(epsilon) || !std::isfinite(r) || r < 0.0) {
            throw std::invalid_argument("Invalid LJ parameters");
        }
        eps_ = epsilon;
        rmin_ = r;
    }

    // Basic setters
    void set_bynu(int bn) { 
        if (bn <= 0) throw std::invalid_argument("Invalid atom number");
        bynu_ = bn; 
    }
    
    void set_type(const std::string& t) { 
        if (!utils::validate::is_valid_name(t, constants::MAX_ATOM_NAME_LENGTH)) {
            throw std::invalid_argument("Invalid atom type");
        }
        type_ = t; 
    }
    
    void set_resname(const std::string& rn) { 
        if (!utils::validate::is_valid_name(rn, constants::MAX_RESIDUE_NAME_LENGTH)) {
            throw std::invalid_argument("Invalid residue name");
        }
        resname_ = rn; 
    }
    
    void set_ires(int ir) { 
        if (ir <= 0) throw std::invalid_argument("Invalid residue number");
        ires_ = ir; 
    }
    
    void set_segid(const std::string& sid) { 
        if (sid.length() > constants::MAX_SEGMENT_NAME_LENGTH) {
            throw std::invalid_argument("Segment ID too long");
        }
        segid_ = sid; 
    }

    // PDB setters
    void set_chain(char ch) { chain_ = ch; }
    void set_hetatm(bool het) { hetatm_ = het; }
    void set_altloc(char alt) { altloc_ = alt; }
    void set_inscode(char ins) { inscode_ = ins; }
    void set_element(const std::string& elem) { element_ = elem; }

    // Utility methods
    bool has_lj_params() const {
        return std::isfinite(eps_) && std::isfinite(rmin_);
    }

    // IValidatable interface
    bool is_valid() const override {
        return bynu_ > 0 && 
               utils::validate::is_valid_name(type_) && 
               utils::validate::is_valid_name(resname_) &&
               ires_ > 0 && 
               utils::validate::is_valid_mass(mass_) && 
               utils::validate::is_valid_charge(charge_) &&
               std::all_of(coor_.begin(), coor_.end(), utils::validate::is_valid_coordinate) &&
               occupancy_ >= 0.0 && occupancy_ <= 1.0 &&
               std::isfinite(tempfactor_);
    }

    // IIdentifiable interface
    std::string get_id() const override {
        return utils::id::generate_atom_id(type_, ires_, segid_);
    }

    void set_id(const std::string& id) override {
        // Parse ID format: segid:resnum:type
        size_t pos1 = id.find(':');
        size_t pos2 = id.find(':', pos1 + 1);
        if (pos1 != std::string::npos && pos2 != std::string::npos) {
            segid_ = id.substr(0, pos1);
            ires_ = std::stoi(id.substr(pos1 + 1, pos2 - pos1 - 1));
            type_ = id.substr(pos2 + 1);
        }
    }

protected:
    // CHARMM standard fields
    int bynu_;                          ///< Atom number (BYNU)
    std::string type_;                  ///< Atom type (TYPE)
    std::string resname_;               ///< Residue name (RESNAME)
    int ires_;                          ///< Residue number (IRES)
    std::string segid_;                 ///< Segment ID (SEGID)
    int iseg_;                          ///< Segment number (ISEG)
    int igro_;                          ///< Group number (IGRO)
    char altloc_;                       ///< Alternate location (PDB)
    char chain_;                        ///< Chain identifier (PDB)
    char inscode_;                      ///< Insertion code (PDB)
    std::array<double, 3> coor_;        ///< Coordinates (X, Y, Z)
    std::array<double, 3> xcom_;        ///< Component coordinates
    std::array<double, 3> xref_;        ///< Reference coordinates
    double occupancy_;                  ///< PDB occupancy
    double tempfactor_;                 ///< PDB temperature factor
    double wmain_;                      ///< Main weight (WMAIN)
    double wcomp_;                      ///< Component weight (WCOMP)
    double mass_;                       ///< Mass (MASS)
    double charge_;                     ///< Charge (CHARGE)
    std::string chem_;                  ///< Chemical type (CHEM)
    double radius_;                     ///< VDW radius (RADIUS)
    double alpha_;                      ///< Polarizability (ALPHA)
    double eps_;                        ///< LJ well depth (epsilon)
    double rmin_;                       ///< LJ Rmin/2 (rmin)
    double fbeta_;                      ///< Force beta (FBETA)
    int move_;                          ///< Movement flag (MOVE)
    int ignore_;                        ///< Ignore flag (IGNORE)
    int constrain_;                     ///< Constraint flag (CONSTRAIN)
    bool hetatm_;                       ///< HETATM flag
    bool initial_;                      ///< Has initial coordinates

    // Scalar properties (SCA1-9)
    std::array<double, 9> scalar_;      ///< User-defined scalar properties

    // Additional PDB fields
    std::string element_;               ///< Element symbol
    std::string chargestr_;             ///< Charge as string from PDB
};

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_ATOM_CORE_HPP 