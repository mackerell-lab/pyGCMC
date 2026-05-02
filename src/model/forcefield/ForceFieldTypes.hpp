#pragma once

#ifndef PYGCMC_MODEL_FORCEFIELD_TYPES_HPP
#define PYGCMC_MODEL_FORCEFIELD_TYPES_HPP

#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <set>

namespace pygcmc {
namespace model {
namespace forcefield {

/**
 * @brief Parameters for non-bonded interactions
 */
struct NonbondedParams {
    int nbxmod = 5;
    bool cdiel = false;
    bool fshift = false;
    bool vatom = false;
    bool vdistance = false;
    bool vfswitch = false;
    double cutnb = 14.0;
    double ctofnb = 12.0;
    double ctonnb = 10.0;
    double eps = 1.0;
    double e14fac = 1.0;
    double wmin = 1.5;
};

/**
 * @brief Parameters for Lennard-Jones interactions
 */
struct LJParams {
    double epsilon = 0.0;     ///< Well depth (kcal/mole)
    double rmin_half = 0.0;   ///< Rmin/2: HALF of the distance at minimum energy (Angstroms)
};

/**
 * @brief Parameters for bond interactions
 */
struct BondParams {
    double kb = 0.0;    ///< Force constant
    double b0 = 0.0;    ///< Equilibrium length
};

/**
 * @brief Parameters for angle interactions
 */
struct AngleParams {
    double ktheta = 0.0;  ///< Force constant
    double theta0 = 0.0;  ///< Equilibrium angle
    double kub = 0.0;     ///< Urey-Bradley force constant
    double s0 = 0.0;      ///< Urey-Bradley equilibrium distance
};

/**
 * @brief Parameters for dihedral interactions
 */
struct DihedralParams {
    double kchi = 0.0;   ///< Force constant
    int n = 1;           ///< Multiplicity
    double delta = 0.0;  ///< Phase shift
};

/**
 * @brief Parameters for improper dihedral interactions
 */
struct ImproperParams {
    double kpsi = 0.0;   ///< Force constant
    double psi0 = 0.0;   ///< Equilibrium angle
};

/**
 * @brief Parameters for NBFIX (specific nonbonded interaction parameters)
 */
struct NBFIXParams {
    double epsilon = 0.0;     ///< Well depth (kcal/mole)
    double rmin = 0.0;        ///< Distance at minimum energy (Angstroms)
};

/**
 * @brief Parameters for Drude polarizability (ALPHA/THOLE)
 */
struct AlphaTHoleParams {
    double alpha = 0.0;       ///< Polarizability (Angstrom^3)
    double thole = 0.0;       ///< Thole screening parameter (dimensionless)
};

/**
 * @brief Parameters for lone pair particles
 */
struct LonePairParams {
    std::string type;         ///< Type of lone pair (e.g., "bisector", "relative")
    std::string host;         ///< Host atom type
    std::string atom1;        ///< First reference atom
    std::string atom2;        ///< Second reference atom
    std::string atom3;        ///< Third reference atom (if needed)
    double distance = 0.0;    ///< Distance from host atom (Angstroms)
    double angle = 0.0;       ///< Angle parameter (degrees)
    double dihedral = 0.0;    ///< Dihedral parameter (degrees)
};

/**
 * @brief Parameters for anisotropic polarizability
 */
struct AnisotropyParams {
    std::string type;         ///< Atom type
    double a11 = 0.0;         ///< XX component of polarizability tensor
    double a22 = 0.0;         ///< YY component of polarizability tensor
    double a33 = 0.0;         ///< ZZ component of polarizability tensor (optional)
    // Note: For 2D anisotropy in Drude model, often only a11 and a22 are used
};

/**
 * @brief Parameters for NBTHOLE (pairwise Thole screening parameters)
 */
struct NBTHOLEParams {
    double thole = 0.0;       ///< Pairwise Thole screening parameter (dimensionless)
};

/**
 * @brief Global Drude model parameters
 */
struct DrudeGlobalParams {
    double tcut = 5.0;         ///< Thole screening cutoff distance (Angstroms)
    int maxnbthole = 5000;     ///< Maximum number of NBTHOLE entries allowed
};



} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_TYPES_HPP
