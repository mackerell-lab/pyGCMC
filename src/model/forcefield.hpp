// src/model/forcefield.hpp

#pragma once

#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <algorithm>
#include <optional>
#include <memory>
#include <set>

namespace pygcmc {
namespace model {

// Forward declarations
class ForceField;

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
 * @brief Detailed explanation of nonbonded parameters in CHARMM force field
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
 * @brief Main force field class that holds all force field parameters
 */
class ForceField {
public:
    ForceField() = default;
    ~ForceField() = default;

    // Add methods
    void add_atom_mass(const std::string& type, double mass) {
        atom_masses_[type] = mass;
    }

    /**
     * @brief Add Lennard-Jones parameters for an atom type
     * @param type The atom type
     * @param epsilon Well depth (kcal/mole)
     * @param rmin_half Rmin/2: HALF of the distance at minimum energy (Angstroms)
     */
    void add_lj_params(const std::string& type, double epsilon, double rmin_half) {
        LJParams params{epsilon, rmin_half};
        lj_params_[type] = params;
    }

    void add_nbfix(const std::string& type1, const std::string& type2, double epsilon) {
        auto key = makeTypePair(type1, type2);
        nbfix_[key] = epsilon;
    }

    void add_bond_params(const std::string& type1, const std::string& type2, double kb, double b0) {
        auto key = makeTypePair(type1, type2);
        BondParams params{kb, b0};
        bond_params_[key] = params;
    }

    void add_angle_params(const std::string& type1, const std::string& type2, const std::string& type3,
                         double ktheta, double theta0, double kub = 0.0, double s0 = 0.0) {
        auto key = makeTypeTriple(type1, type2, type3);
        AngleParams params{ktheta, theta0, kub, s0};
        angle_params_[key] = params;
    }

    void add_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kchi, int n, double delta) {
        auto key = makeTypeQuad(type1, type2, type3, type4);
        DihedralParams params{kchi, n, delta};
        dihedral_params_[key].push_back(params);
    }

    void add_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4,
                            double kpsi, double psi0) {
        auto key = makeTypeQuad(type1, type2, type3, type4);
        ImproperParams params{kpsi, psi0};
        improper_params_[key] = params;
    }

    // Getter methods with error checking
    double get_atom_mass(const std::string& type) const {
        auto it = atom_masses_.find(type);
        if (it == atom_masses_.end()) {
            throw std::runtime_error("Atom mass not found for type: " + type);
        }
        return it->second;
    }

    const LJParams& get_lj_params(const std::string& type) const {
        auto it = lj_params_.find(type);
        if (it == lj_params_.end()) {
            throw std::runtime_error("LJ parameters not found for type: " + type);
        }
        return it->second;
    }

    std::pair<double, bool> get_nbfix(const std::string& type1, const std::string& type2) const {
        auto key = makeTypePair(type1, type2);
        auto it = nbfix_.find(key);
        if (it == nbfix_.end()) {
            return std::make_pair(0.0, false);
        }
        return std::make_pair(it->second, true);
    }

    const BondParams& get_bond_params(const std::string& type1, const std::string& type2) const {
        auto key = makeTypePair(type1, type2);
        auto it = bond_params_.find(key);
        if (it == bond_params_.end()) {
            throw std::runtime_error("Bond parameters not found for types: " + type1 + "-" + type2);
        }
        return it->second;
    }

    const AngleParams& get_angle_params(const std::string& type1,
                                      const std::string& type2,
                                      const std::string& type3) const {
        // Try both orientations of the outer atoms while keeping the middle atom fixed
        auto key1 = std::make_tuple(type1, type2, type3);
        auto it = angle_params_.find(key1);
        if (it != angle_params_.end()) {
            return it->second;
        }

        // Try the reverse orientation
        auto key2 = std::make_tuple(type3, type2, type1);
        it = angle_params_.find(key2);
        if (it != angle_params_.end()) {
            return it->second;
        }

        throw std::runtime_error("Angle parameters not found for types: " + 
                               type1 + "-" + type2 + "-" + type3);
    }

    const std::vector<DihedralParams>& get_dihedral_params(const std::string& type1,
                                                          const std::string& type2,
                                                          const std::string& type3,
                                                          const std::string& type4) const {
        auto key = makeTypeQuad(type1, type2, type3, type4);
        auto it = dihedral_params_.find(key);
        if (it == dihedral_params_.end()) {
            throw std::runtime_error("Dihedral parameters not found for types: " +
                                   type1 + "-" + type2 + "-" + type3 + "-" + type4);
        }
        return it->second;
    }

    const ImproperParams& get_improper_params(const std::string& type1, const std::string& type2,
                                            const std::string& type3, const std::string& type4) const {
        auto key = makeTypeQuad(type1, type2, type3, type4);
        auto it = improper_params_.find(key);
        if (it == improper_params_.end()) {
            throw std::runtime_error("Improper parameters not found for types: " +
                                   type1 + "-" + type2 + "-" + type3 + "-" + type4);
        }
        return it->second;
    }

    // Existence check methods
    bool has_atom_mass(const std::string& type) const {
        return atom_masses_.find(type) != atom_masses_.end();
    }

    bool has_lj_params(const std::string& type) const {
        return lj_params_.find(type) != lj_params_.end();
    }

    bool has_nbfix(const std::string& type1, const std::string& type2) const {
        auto key = makeTypePair(type1, type2);
        return nbfix_.find(key) != nbfix_.end();
    }

    bool has_bond_params(const std::string& type1, const std::string& type2) const {
        auto key = makeTypePair(type1, type2);
        return bond_params_.find(key) != bond_params_.end();
    }

    bool has_angle_params(const std::string& type1, const std::string& type2,
                         const std::string& type3) const {
        auto key = makeTypeTriple(type1, type2, type3);
        return angle_params_.find(key) != angle_params_.end();
    }

    bool has_dihedral_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const {
        auto key = makeTypeQuad(type1, type2, type3, type4);
        return dihedral_params_.find(key) != dihedral_params_.end();
    }

    bool has_improper_params(const std::string& type1, const std::string& type2,
                            const std::string& type3, const std::string& type4) const {
        auto key = makeTypeQuad(type1, type2, type3, type4);
        return improper_params_.find(key) != improper_params_.end();
    }

    // Size methods
    size_t get_num_atom_types() const { return atom_masses_.size(); }
    size_t get_num_lj_params() const { return lj_params_.size(); }
    size_t get_num_nbfix() const { return nbfix_.size(); }
    size_t get_num_bond_types() const { return bond_params_.size(); }
    size_t get_num_angle_types() const { return angle_params_.size(); }
    size_t get_num_dihedral_types() const { return dihedral_params_.size(); }
    size_t get_num_improper_types() const { return improper_params_.size(); }

    // Access to nonbonded parameters
    const NonbondedParams& get_nonbonded_params() const { return nonbonded_params_; }
    NonbondedParams& get_nonbonded_params() { return nonbonded_params_; }

    // Static helper methods for making parameter keys
    static std::pair<std::string, std::string> makeTypePair(const std::string& type1,
                                                           const std::string& type2) {
        return type1 < type2 ? std::make_pair(type1, type2) : std::make_pair(type2, type1);
    }

    static std::tuple<std::string, std::string, std::string> makeTypeTriple(
        const std::string& type1, const std::string& type2, const std::string& type3) {
        // For angle parameters in CHARMM force field:
        // 1. The middle atom (type2) must stay in the middle
        // 2. Try both orientations of the outer atoms
        // Return both (type1, type2, type3) and (type3, type2, type1)
        // This ensures we find the parameter regardless of how it's stored in the force field
        return std::make_tuple(type1, type2, type3);
    }

    static std::tuple<std::string, std::string, std::string, std::string> makeTypeQuad(
        const std::string& type1, const std::string& type2,
        const std::string& type3, const std::string& type4) {
        return std::make_tuple(type1, type2, type3, type4);
    }

    // Direct access to parameter maps (for Python bindings)
    const std::map<std::string, double>& get_atom_masses() const { return atom_masses_; }
    const std::map<std::string, LJParams>& get_lj_params() const { return lj_params_; }
    const std::map<std::pair<std::string, std::string>, double>& get_nbfix() const { return nbfix_; }
    const std::map<std::pair<std::string, std::string>, BondParams>& get_bond_params() const { return bond_params_; }
    const std::map<std::tuple<std::string, std::string, std::string>, AngleParams>& get_angle_params() const { return angle_params_; }
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>>& get_dihedral_params() const { return dihedral_params_; }
    const std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams>& get_improper_params() const { return improper_params_; }

private:
    // Parameter storage
    std::map<std::string, double> atom_masses_;
    std::map<std::string, LJParams> lj_params_;
    std::map<std::pair<std::string, std::string>, double> nbfix_;
    std::map<std::pair<std::string, std::string>, BondParams> bond_params_;
    std::map<std::tuple<std::string, std::string, std::string>, AngleParams> angle_params_;
    std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>> dihedral_params_;
    std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams> improper_params_;
    NonbondedParams nonbonded_params_;
};

} // namespace model
} // namespace pygcmc


