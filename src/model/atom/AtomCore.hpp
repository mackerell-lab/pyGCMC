#pragma once

#ifndef PYGCMC_MODEL_ATOM_CORE_HPP
#define PYGCMC_MODEL_ATOM_CORE_HPP

#include "../common/ModelInterface.hpp"
#include "../common/ModelConstants.hpp"
#include <string>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <algorithm>

namespace pygcmc {
namespace model {
namespace atom {

/**
 * @brief Core Atom class following CHARMM naming conventions and PDB format
 * This is the minimal, core implementation with essential functionality
 */
class AtomCore : public common::IValidatable, public common::IIdentifiable {
public:
    // Core constructors
    AtomCore() : 
        bynu(0), ires(0), iseg(0), igro(0),
        altloc(common::constants::DEFAULT_ALTLOC),
        chain(common::constants::DEFAULT_CHAIN),
        inscode(common::constants::DEFAULT_INSCODE),
        coor{0.0, 0.0, 0.0},
        occupancy(common::constants::DEFAULT_OCCUPANCY),
        tempfactor(common::constants::DEFAULT_TEMPFACTOR),
        wmain(common::constants::DEFAULT_WMAIN),
        wcomp(common::constants::DEFAULT_WCOMP),
        mass(0.0), charge(0.0), radius(0.0), alpha(0.0),
        eps(common::constants::INVALID_DOUBLE),
        rmin(common::constants::INVALID_DOUBLE),
        fbeta(0.0),
        move(common::constants::DEFAULT_MOVE),
        ignore(common::constants::DEFAULT_IGNORE),
        constrain(common::constants::DEFAULT_CONSTRAIN),
        hetatm(false), initial(false)
    {}

    AtomCore(int bynu, const std::string& type, const std::string& resname,
             int ires, const std::string& segid = "", int iseg = 0,
             double x = 0.0, double y = 0.0, double z = 0.0,
             double wmain = common::constants::DEFAULT_WMAIN, 
             double mass = 0.0, double charge = 0.0,
             const std::string& chem = "", bool hetatm = false) :
        bynu(bynu), type(type), resname(resname), ires(ires),
        segid(segid), iseg(iseg), igro(0),
        altloc(common::constants::DEFAULT_ALTLOC),
        chain(common::constants::DEFAULT_CHAIN),
        inscode(common::constants::DEFAULT_INSCODE),
        coor{x, y, z},
        occupancy(common::constants::DEFAULT_OCCUPANCY),
        tempfactor(common::constants::DEFAULT_TEMPFACTOR),
        wmain(wmain), wcomp(common::constants::DEFAULT_WCOMP),
        mass(mass), charge(charge), chem(chem),
        eps(common::constants::INVALID_DOUBLE),
        rmin(common::constants::INVALID_DOUBLE),
        hetatm(hetatm), initial(false) {}

    virtual ~AtomCore() = default;

    // IIdentifiable interface
    int get_id() const override { return bynu; }
    void set_id(int id) override { 
        if (id <= 0) throw std::invalid_argument("Invalid atom number");
        bynu = id; 
    }

    // IValidatable interface
    bool is_valid() const override {
        return bynu > 0 && !type.empty() && !resname.empty() &&
               ires > 0 && std::isfinite(mass) && std::isfinite(charge) &&
               std::all_of(coor.begin(), coor.end(), 
                          [](double x) { return std::isfinite(x); }) &&
               std::isfinite(occupancy) && std::isfinite(tempfactor);
    }

    // Core getters
    int get_bynu() const noexcept { return bynu; }
    const std::string& get_type() const noexcept { return type; }
    const std::string& get_resname() const noexcept { return resname; }
    int get_ires() const noexcept { return ires; }
    const std::string& get_segid() const noexcept { return segid; }
    int get_iseg() const noexcept { return iseg; }
    int get_igro() const noexcept { return igro; }
    const std::string& get_chem() const noexcept { return chem; }
    double get_wmain() const noexcept { return wmain; }
    double get_mass() const noexcept { return mass; }
    double get_charge() const noexcept { return charge; }
    
    // Coordinate access
    const std::array<double, 3>& get_coor() const noexcept { return coor; }
    double get_x() const noexcept { return coor[0]; }
    double get_y() const noexcept { return coor[1]; }
    double get_z() const noexcept { return coor[2]; }

    // Force field parameters
    double get_eps() const noexcept { return eps; }
    double get_rmin() const noexcept { return rmin; }

    // PDB specific getters
    char get_altloc() const noexcept { return altloc; }
    char get_inscode() const noexcept { return inscode; }
    char get_chain() const noexcept { return chain; }
    bool is_hetatm() const noexcept { return hetatm; }
    double get_occupancy() const noexcept { return occupancy; }
    double get_tempfactor() const noexcept { return tempfactor; }

    // Core setters with validation
    void set_coor(double x, double y, double z) {
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
            throw std::invalid_argument("Invalid coordinates");
        }
        coor = {x, y, z};
    }

    void set_mass_charge(double m, double q) {
        if (!std::isfinite(m) || m < 0.0 || !std::isfinite(q)) {
            throw std::invalid_argument("Invalid mass or charge");
        }
        mass = m;
        charge = q;
    }

    void set_lj_params(double epsilon, double r) {
        if (!std::isfinite(epsilon) || !std::isfinite(r) || r < 0.0) {
            throw std::invalid_argument("Invalid LJ parameters");
        }
        eps = epsilon;
        rmin = r;
    }

    // Basic setters
    void set_type(const std::string& t) { 
        if (t.empty()) throw std::invalid_argument("Empty atom type");
        type = t; 
    }
    
    void set_resname(const std::string& rn) { 
        if (rn.empty()) throw std::invalid_argument("Empty residue name");
        resname = rn; 
    }
    
    void set_ires(int ir) { 
        if (ir <= 0) throw std::invalid_argument("Invalid residue number");
        ires = ir; 
    }
    
    void set_segid(const std::string& sid) { segid = sid; }
    void set_chain(char ch) { chain = ch; }
    void set_hetatm(bool het) { hetatm = het; }
    
    void set_bynu(int num) {
        if (num <= 0) throw std::invalid_argument("Invalid atom number");
        bynu = num;
    }

    // Utility methods
    bool has_lj_params() const {
        return std::isfinite(eps) && std::isfinite(rmin);
    }

protected:
    // CHARMM standard fields
    int bynu;                    ///< Atom number (BYNU)
    std::string type;            ///< Atom type (TYPE)
    std::string resname;         ///< Residue name (RESNAME)
    int ires;                    ///< Residue number (IRES)
    std::string segid;           ///< Segment ID (SEGID)
    int iseg;                    ///< Segment number (ISEG)  
    int igro;                    ///< Group number (IGRO)
    char altloc;                 ///< Alternate location (PDB)
    char chain;                  ///< Chain identifier (PDB)
    char inscode;                ///< Insertion code (PDB)
    std::array<double, 3> coor;  ///< Coordinates (X, Y, Z)
    double occupancy;            ///< PDB occupancy
    double tempfactor;           ///< PDB temperature factor
    double wmain;                ///< Main weight (WMAIN)
    double wcomp;                ///< Component weight (WCOMP)
    double mass;                 ///< Mass (MASS)
    double charge;               ///< Charge (CHARGE)
    std::string chem;            ///< Chemical type (CHEM)
    double radius;               ///< VDW radius (RADIUS)
    double alpha;                ///< Polarizability (ALPHA)
    double eps;                  ///< LJ well depth (epsilon)
    double rmin;                 ///< LJ Rmin/2 (rmin)
    double fbeta;                ///< Force beta (FBETA)
    int move;                    ///< Movement flag (MOVE)
    int ignore;                  ///< Ignore flag (IGNORE)
    int constrain;               ///< Constraint flag (CONSTRAIN)
    bool hetatm;                 ///< HETATM flag
    bool initial;                ///< Has initial coordinates

    // Additional arrays
    std::array<double, 3> xcom;  ///< Component coordinates
    std::array<double, 3> xref;  ///< Reference coordinates
    std::array<double, 9> scalar; ///< User-defined scalar properties
    std::string element;         ///< Element symbol
    std::string chargestr;       ///< Charge as string from PDB
};

} // namespace atom
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_ATOM_CORE_HPP 