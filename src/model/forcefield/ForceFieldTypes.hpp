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
 * 
 * Detailed explanation of nonbonded parameters in CHARMM force field:
 * 
 * The nonbonded energy function in CHARMM consists of Lennard-Jones (LJ) and electrostatic terms:
 * 
 * V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
 * where:
 * - epsilon (Eps,i,j) = sqrt(eps,i * eps,j) [kcal/mole]
 * - Rmin,i,j = Rmin/2,i + Rmin/2,j [Angstroms]
 * - ri,j is the distance between atoms i and j
 * 
 * Parameters:
 * @param nbxmod   Nonbonded exclusion model (5 = use switching functions)
 * @param cdiel    Use constant dielectric (true/false)
 * @param fshift   Use force shifting (true/false)
 * @param vatom    Use atom-based potential (true/false)
 * @param vdistance Use distance-based potential (true/false)
 * @param vfswitch Use force switching (true/false)
 * @param cutnb    Nonbonded cutoff distance [Angstroms]
 * @param ctofnb   Distance at which switching function takes effect for nonbonded [Angstroms]
 * @param ctonnb   Distance at which switching function takes effect for 1-4 interactions [Angstroms]
 * @param eps      Dielectric constant
 * @param e14fac   Scaling factor for 1-4 interactions
 * @param wmin     Minimum weighting in switching function
 * 
 * Units:
 * - Distances in Angstroms
 * - Energies in kcal/mole
 * - Dielectric constant is dimensionless
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
 * 
 * The Lennard-Jones potential is defined as:
 * V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
 * where:
 * - Eps,i,j = sqrt(eps,i * eps,j)
 * - Rmin,i,j = Rmin/2,i + Rmin/2,j
 * 
 * Units:
 * - epsilon: kcal/mole
 * - rmin_half: Angstroms (Rmin/2: HALF of the distance at minimum energy)
 */
struct LJParams {
    double epsilon = 0.0;     ///< Well depth (kcal/mole)
    double rmin_half = 0.0;   ///< Rmin/2: HALF of the distance at minimum energy (Angstroms)
    
    LJParams() = default;
    LJParams(double eps, double rmin_h) : epsilon(eps), rmin_half(rmin_h) {}
};

/**
 * @brief Parameters for bond interactions
 */
struct BondParams {
    double kb = 0.0;    ///< Force constant
    double b0 = 0.0;    ///< Equilibrium length
    
    BondParams() = default;
    BondParams(double kb_val, double b0_val) : kb(kb_val), b0(b0_val) {}
};

/**
 * @brief Parameters for angle interactions
 */
struct AngleParams {
    double ktheta = 0.0;  ///< Force constant
    double theta0 = 0.0;  ///< Equilibrium angle
    double kub = 0.0;     ///< Urey-Bradley force constant
    double s0 = 0.0;      ///< Urey-Bradley equilibrium distance
    
    AngleParams() = default;
    AngleParams(double ktheta_val, double theta0_val, double kub_val = 0.0, double s0_val = 0.0) 
        : ktheta(ktheta_val), theta0(theta0_val), kub(kub_val), s0(s0_val) {}
};

/**
 * @brief Parameters for dihedral interactions
 */
struct DihedralParams {
    double kchi = 0.0;   ///< Force constant
    int n = 1;           ///< Multiplicity
    double delta = 0.0;  ///< Phase shift
    
    DihedralParams() = default;
    DihedralParams(double kchi_val, int n_val, double delta_val) 
        : kchi(kchi_val), n(n_val), delta(delta_val) {}
};

/**
 * @brief Parameters for improper dihedral interactions
 */
struct ImproperParams {
    double kpsi = 0.0;   ///< Force constant
    double psi0 = 0.0;   ///< Equilibrium angle
    
    ImproperParams() = default;
    ImproperParams(double kpsi_val, double psi0_val) : kpsi(kpsi_val), psi0(psi0_val) {}
};

/**
 * @brief Parameters for NBFIX (specific nonbonded interaction parameters)
 * 
 * In CHARMM, NBFIX allows specification of specific Lennard-Jones parameters
 * for particular pairs of atom types, overriding the standard combining rules.
 * 
 * Units:
 * - epsilon: kcal/mole
 * - rmin: Angstroms (full Rmin, not Rmin/2)
 */
struct NBFIXParams {
    double epsilon = 0.0;     ///< Well depth (kcal/mole)
    double rmin = 0.0;        ///< Distance at minimum energy (Angstroms)
    
    NBFIXParams() = default;
    NBFIXParams(double eps, double rmin_val) : epsilon(eps), rmin(rmin_val) {}
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
    /**
     * @brief Create a parameter key for pair interactions
     * Orders the types alphabetically for consistent lookup
     */
    static std::pair<std::string, std::string> makeTypePair(const std::string& type1,
                                                           const std::string& type2) {
        return type1 < type2 ? std::make_pair(type1, type2) : std::make_pair(type2, type1);
    }

    /**
     * @brief Create a parameter key for angle interactions
     * For angle parameters in CHARMM force field:
     * 1. The middle atom (type2) must stay in the middle
     * 2. Store parameters in the order they appear in the parameter file
     */
    static std::tuple<std::string, std::string, std::string> makeTypeTriple(
        const std::string& type1, const std::string& type2, const std::string& type3) {
        return std::make_tuple(type1, type2, type3);
    }

    /**
     * @brief Create a parameter key for dihedral/improper interactions
     */
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