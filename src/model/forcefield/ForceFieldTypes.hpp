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
 * @brief Statistics for force field parameters
 */
struct ForceFieldStats {
    size_t num_atom_types = 0;
    size_t num_lj_params = 0;
    size_t num_nbfix = 0;
    size_t num_bond_types = 0;
    size_t num_angle_types = 0;
    size_t num_dihedral_types = 0;
    size_t num_improper_types = 0;
};

/**
 * @brief Completeness check result
 */
struct CompletenessResult {
    bool is_complete = false;
    std::set<std::string> missing_atom_masses;
    std::set<std::string> missing_lj_params;
    std::string summary;
};

/**
 * @brief Utility class for creating parameter keys
 */
class ParamKeyUtils {
public:
    static std::pair<std::string, std::string> makeTypePair(const std::string& type1, const std::string& type2) {
        return type1 < type2 ? std::make_pair(type1, type2) : std::make_pair(type2, type1);
    }

    static std::tuple<std::string, std::string, std::string> makeTypeTriple(
        const std::string& type1, const std::string& type2, const std::string& type3) {
        return std::make_tuple(type1, type2, type3);
    }

    static std::tuple<std::string, std::string, std::string, std::string> makeTypeQuad(
        const std::string& type1, const std::string& type2,
        const std::string& type3, const std::string& type4) {
        return std::make_tuple(type1, type2, type3, type4);
    }
};

} // namespace forcefield
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_FORCEFIELD_TYPES_HPP 