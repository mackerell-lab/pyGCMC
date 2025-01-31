// src/model/atom.hpp

#pragma once

#ifndef PYGCMC_MODEL_ATOM_HPP
#define PYGCMC_MODEL_ATOM_HPP

#include <string>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <algorithm>  // for std::all_of

namespace pygcmc {
namespace model {

/**
 * @brief Atom class following CHARMM naming conventions and PDB format
 */
class Atom {
public:
    // Constructors
    Atom() : 
        bynu(0),           // BYNU - atom number
        ires(0),           // IRES - residue number
        iseg(0),           // ISEG - segment number
        igro(0),           // IGRO - group number
        altloc(' '),       // Alternate location (PDB specific)
        chain(' '),        // Chain identifier (PDB specific)
        inscode(' '),      // Insertion code (PDB specific)
        coor{0.0, 0.0, 0.0}, // Coordinates
        occupancy(1.0),    // PDB occupancy
        tempfactor(0.0),   // PDB temperature factor
        wmain(1.0),        // Main weight (WMAIN)
        wcomp(1.0),        // Component weight (WCOMP)
        mass(0.0),         // Mass (MASS)
        charge(0.0),       // Charge (CHARGE)
        radius(0.0),       // VDW radius (RADIUS)
        alpha(0.0),        // Polarizability (ALPHA)
        eps(std::numeric_limits<double>::quiet_NaN()),  // LJ epsilon
        rmin(std::numeric_limits<double>::quiet_NaN()), // LJ Rmin/2
        fbeta(0.0),        // Force beta (FBETA)
        move(1),           // Movement flag (MOVE)
        ignore(0),         // Ignore flag (IGNORE)
        constrain(0),      // Constraint flag (CONSTRAIN)
        hetatm(false),     // HETATM flag
        initial(false)     // Has initial coordinates
    {}

    Atom(int bynu, const std::string& type, const std::string& resname,
         int ires, const std::string& segid = "", int iseg = 0,
         double x = 0.0, double y = 0.0, double z = 0.0,
         double wmain = 1.0, double mass = 0.0, double charge = 0.0,
         const std::string& chem = "", bool hetatm = false) :
        bynu(bynu),
        type(type),
        resname(resname),
        ires(ires),
        segid(segid),
        iseg(iseg),
        igro(0),
        altloc(' '),
        chain(' '),
        inscode(' '),
        coor{x, y, z},
        wmain(wmain),
        mass(mass),
        charge(charge),
        chem(chem),
        eps(std::numeric_limits<double>::quiet_NaN()),
        rmin(std::numeric_limits<double>::quiet_NaN()),
        hetatm(hetatm) {}

    // CHARMM standard getters
    int get_bynu() const noexcept { return bynu; }              // Atom number
    const std::string& get_type() const noexcept { return type; }    // Atom type
    const std::string& get_resname() const noexcept { return resname; } // Residue name
    int get_ires() const noexcept { return ires; }              // Residue number
    const std::string& get_segid() const noexcept { return segid; }   // Segment ID
    int get_iseg() const noexcept { return iseg; }              // Segment number
    int get_igro() const noexcept { return igro; }              // Group number
    const std::string& get_chem() const noexcept { return chem; }     // Chemical type
    double get_wmain() const noexcept { return wmain; }         // Main weight
    double get_mass() const noexcept { return mass; }           // Mass
    double get_charge() const noexcept { return charge; }       // Charge
    
    // Coordinate access
    const std::array<double, 3>& get_coor() const noexcept { return coor; }
    double get_x() const noexcept { return coor[0]; }
    double get_y() const noexcept { return coor[1]; }
    double get_z() const noexcept { return coor[2]; }

    // Force field parameters
    double get_eps() const noexcept { return eps; }    // LJ well depth
    double get_rmin() const noexcept { return rmin; }  // LJ Rmin/2

    // PDB specific getters
    char get_altloc() const noexcept { return altloc; }
    char get_inscode() const noexcept { return inscode; }
    char get_chain() const noexcept { return chain; }
    bool is_hetatm() const noexcept { return hetatm; }

    // Additional getters for PDB compatibility
    double get_occupancy() const noexcept { return occupancy; }
    double get_tempfactor() const noexcept { return tempfactor; }
    const std::string& get_element() const noexcept { return element; }
    const std::string& get_charge_string() const noexcept { return chargestr; }

    // Setters with validation
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

    // Additional PDB setters
    void set_chain(char ch) { chain = ch; }
    void set_hetatm(bool het) { hetatm = het; }

    // Additional setters
    void set_bynu(int bn) { 
        if (bn <= 0) throw std::invalid_argument("Invalid atom number");
        bynu = bn; 
    }
    
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
    
    void set_altloc(char alt) { altloc = alt; }
    void set_inscode(char ins) { inscode = ins; }

    // Utility methods
    bool has_lj_params() const {
        return std::isfinite(eps) && std::isfinite(rmin);
    }

    bool is_valid() const {
        return bynu > 0 && !type.empty() && !resname.empty() &&
               ires > 0 && std::isfinite(mass) && std::isfinite(charge) &&
               std::all_of(coor.begin(), coor.end(), 
                          [](double x) { return std::isfinite(x); }) &&
               std::isfinite(occupancy) && std::isfinite(tempfactor);
    }

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

private:
    // CHARMM standard fields
    int bynu;                  ///< Atom number (BYNU)
    std::string type;          ///< Atom type (TYPE)
    std::string resname;       ///< Residue name (RESNAME)
    int ires;                  ///< Residue number (IRES)
    std::string segid;         ///< Segment ID (SEGID)
    int iseg;                  ///< Segment number (ISEG)
    int igro;                  ///< Group number (IGRO)
    char altloc;               ///< Alternate location (PDB)
    char chain;                ///< Chain identifier (PDB)
    char inscode;              ///< Insertion code (PDB)
    std::array<double, 3> coor; ///< Coordinates (X, Y, Z)
    std::array<double, 3> xcom; ///< Component coordinates (XCOM, YCOM, ZCOM)
    std::array<double, 3> xref; ///< Reference coordinates (XREF, YREF, ZREF)
    double occupancy;          ///< PDB occupancy
    double tempfactor;         ///< PDB temperature factor
    double wmain;              ///< Main weight (WMAIN)
    double wcomp;              ///< Component weight (WCOMP)
    double mass;               ///< Mass (MASS)
    double charge;             ///< Charge (CHARGE)
    std::string chem;          ///< Chemical type (CHEM)
    double radius;             ///< VDW radius (RADIUS)
    double alpha;              ///< Polarizability (ALPHA)
    double eps;                ///< LJ well depth (epsilon)
    double rmin;               ///< LJ Rmin/2 (rmin)
    double fbeta;              ///< Force beta (FBETA)
    int move;                  ///< Movement flag (MOVE)
    int ignore;               ///< Ignore flag (IGNORE)
    int constrain;            ///< Constraint flag (CONSTRAIN)
    bool hetatm;              ///< HETATM flag
    bool initial;             ///< Has initial coordinates

    // Scalar properties (SCA1-9)
    std::array<double, 9> scalar;  ///< User-defined scalar properties

    // Additional PDB fields
    std::string element;       ///< Element symbol
    std::string chargestr;     ///< Charge as string from PDB
};

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_ATOM_HPP

