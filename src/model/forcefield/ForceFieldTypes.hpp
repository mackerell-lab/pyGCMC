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



} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_TYPES_HPP 