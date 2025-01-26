// src/data/atom.hpp

#ifndef PYGCMC_DATA_ATOM_HPP
#define PYGCMC_DATA_ATOM_HPP

#include <string>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <algorithm>  // for std::all_of

namespace pygcmc {
namespace data {

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
    int getBynu() const noexcept { return bynu; }              // Atom number
    const std::string& getType() const noexcept { return type; }    // Atom type
    const std::string& getResname() const noexcept { return resname; } // Residue name
    int getIres() const noexcept { return ires; }              // Residue number
    const std::string& getSegid() const noexcept { return segid; }   // Segment ID
    int getIseg() const noexcept { return iseg; }              // Segment number
    int getIgro() const noexcept { return igro; }              // Group number
    const std::string& getChem() const noexcept { return chem; }     // Chemical type
    double getWmain() const noexcept { return wmain; }         // Main weight
    double getMass() const noexcept { return mass; }           // Mass
    double getCharge() const noexcept { return charge; }       // Charge
    
    // Coordinate access
    const std::array<double, 3>& getCoor() const noexcept { return coor; }
    double getX() const noexcept { return coor[0]; }
    double getY() const noexcept { return coor[1]; }
    double getZ() const noexcept { return coor[2]; }

    // Force field parameters
    double getEps() const noexcept { return eps; }    // LJ well depth
    double getRmin() const noexcept { return rmin; }  // LJ Rmin/2

    // PDB specific getters
    char getAltloc() const noexcept { return altloc; }
    char getInscode() const noexcept { return inscode; }
    char getChain() const noexcept { return chain; }
    bool isHetatm() const noexcept { return hetatm; }

    // Additional getters for PDB compatibility
    double getOccupancy() const noexcept { return occupancy; }
    double getTempfactor() const noexcept { return tempfactor; }
    const std::string& getElement() const noexcept { return element; }
    const std::string& getChargeString() const noexcept { return chargestr; }

    // Setters with validation
    void setCoor(double x, double y, double z) {
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
            throw std::invalid_argument("Invalid coordinates");
        }
        coor = {x, y, z};
    }

    void setMassCharge(double m, double q) {
        if (!std::isfinite(m) || m < 0.0 || !std::isfinite(q)) {
            throw std::invalid_argument("Invalid mass or charge");
        }
        mass = m;
        charge = q;
    }

    void setLJParams(double epsilon, double r) {
        if (!std::isfinite(epsilon) || !std::isfinite(r) || r < 0.0) {
            throw std::invalid_argument("Invalid LJ parameters");
        }
        eps = epsilon;
        rmin = r;
    }

    // Additional setters for PDB compatibility
    void setOccupancy(double occ) {
        if (occ < 0.0 || occ > 1.0) {
            throw std::invalid_argument("Occupancy must be between 0 and 1");
        }
        occupancy = occ;
    }

    void setTempfactor(double temp) {
        if (!std::isfinite(temp)) {
            throw std::invalid_argument("Invalid temperature factor");
        }
        tempfactor = temp;
    }

    void setElement(const std::string& elem) {
        element = elem;
    }

    void setChargeString(const std::string& chg) {
        chargestr = chg;
    }

    // Additional PDB setters
    void setChain(char ch) { chain = ch; }
    void setHetatm(bool het) { hetatm = het; }

    // Utility methods
    bool hasLJParams() const {
        return std::isfinite(eps) && std::isfinite(rmin);
    }

    bool isValid() const {
        return bynu > 0 && !type.empty() && !resname.empty() &&
               ires > 0 && std::isfinite(mass) && std::isfinite(charge) &&
               std::all_of(coor.begin(), coor.end(), 
                          [](double x) { return std::isfinite(x); }) &&
               std::isfinite(occupancy) && std::isfinite(tempfactor);
    }

    // PDB format utilities
    static std::string formatPDBAtomName(const std::string& name) {
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

    std::string getFormattedAtomName() const {
        return formatPDBAtomName(type);
    }

    std::string getResidueID() const {
        // Combine residue number and insertion code (e.g., "153A")
        if (inscode == ' ') {
            return std::to_string(ires);
        }
        return std::to_string(ires) + inscode;
    }

    void setResidueID(const std::string& resid) {
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
    double wmain;              ///< Main weight (WMAIN)
    double wcomp;              ///< Component weight (WCOMP)
    double mass;               ///< Mass (MASS)
    double charge;             ///< Charge (CHARGE)
    std::string chem;          ///< Chemical type (CHEM)
    double eps;                ///< LJ well depth (epsilon)
    double rmin;               ///< LJ Rmin/2 (rmin)
    double radius;             ///< VDW radius (RADIUS)
    double alpha;              ///< Polarizability (ALPHA)
    double fbeta;              ///< Force beta (FBETA)
    int move;                  ///< Movement flag (MOVE)
    int ignore;               ///< Ignore flag (IGNORE)
    int constrain;            ///< Constraint flag (CONSTRAIN)
    bool hetatm;              ///< HETATM flag
    bool initial;             ///< Has initial coordinates

    // Scalar properties (SCA1-9)
    std::array<double, 9> scalar;  ///< User-defined scalar properties

    // Additional PDB fields
    double occupancy;          ///< PDB occupancy
    double tempfactor;         ///< PDB temperature factor
    std::string element;       ///< Element symbol
    std::string chargestr;     ///< Charge as string from PDB
};

} // namespace data
} // namespace pygcmc

#endif // PYGCMC_DATA_ATOM_HPP

