// src/model/forcefield.hpp

#pragma once

#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <stdexcept>
#include <algorithm>  // for std::min_element

namespace pygcmc {

// 非键相互作用参数
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

// LJ参数
struct LJParams {
    double epsilon = 0.0;
    double rmin = 0.0;
};

// 键参数
struct BondParams {
    double kb = 0.0;    // force constant
    double b0 = 0.0;    // equilibrium length
};

// 角度参数
struct AngleParams {
    double ktheta = 0.0;  // force constant
    double theta0 = 0.0;  // equilibrium angle
    double kub = 0.0;     // Urey-Bradley force constant
    double s0 = 0.0;      // Urey-Bradley equilibrium distance
};

// 二面角参数
struct DihedralParams {
    double kchi = 0.0;   // force constant
    int n = 1;           // multiplicity
    double delta = 0.0;  // phase shift
};

// 非正常二面角参数
struct ImproperParams {
    double kpsi = 0.0;   // force constant
    double psi0 = 0.0;   // equilibrium angle
};

class ForceField {
public:
    // 存储结构
    std::map<std::string, double> atom_masses;                    // 原子质量
    std::map<std::string, LJParams> lj_params;                    // LJ参数
    std::map<std::pair<std::string, std::string>, double> nbfix;  // NBFIX参数
    std::map<std::pair<std::string, std::string>, BondParams> bond_params;  // 键参数
    std::map<std::tuple<std::string, std::string, std::string>, AngleParams> angle_params;  // 角度参数
    std::map<std::tuple<std::string, std::string, std::string, std::string>, 
            std::vector<DihedralParams>> dihedral_params;  // 二面角参数
    std::map<std::tuple<std::string, std::string, std::string, std::string>, 
            ImproperParams> improper_params;  // 非正常二面角参数
    NonbondedParams nonbonded_params;  // 非键相互作用参数

    // Helper functions for making keys
    static std::pair<std::string, std::string> makeTypePair(const std::string& type1, const std::string& type2) {
        return type1 < type2 ? std::make_pair(type1, type2) : std::make_pair(type2, type1);
    }

    static std::tuple<std::string, std::string, std::string> makeTypeTriple(
        const std::string& type1, const std::string& type2, const std::string& type3) {
        // For angle parameters, we need to handle both symmetric and alternative representations
        // Create a vector of all possible representations
        std::vector<std::tuple<std::string, std::string, std::string>> keys = {
            std::make_tuple(type1, type2, type3),  // original order
            std::make_tuple(type3, type2, type1),  // symmetric order
            std::make_tuple(type2, type1, type3),  // alternative representation
            std::make_tuple(type2, type3, type1)   // symmetric alternative representation
        };
        
        // Return the lexicographically smallest key to ensure consistency
        return *std::min_element(keys.begin(), keys.end());
    }

    static std::tuple<std::string, std::string, std::string, std::string> makeTypeQuad(
        const std::string& type1, const std::string& type2, const std::string& type3, const std::string& type4) {
        // For dihedral parameters in CHARMM, we need to preserve the order as defined in the parameter file
        // The order in the parameter file is the correct one, we should not change it
        return std::make_tuple(type1, type2, type3, type4);
    }

    // Getter methods
    const NonbondedParams& getNonbondedParams() const {
        return nonbonded_params;
    }

    const LJParams& getLJParams(const std::string& type) const {
        auto it = lj_params.find(type);
        if (it == lj_params.end()) {
            throw std::runtime_error("LJ parameters not found for type: " + type);
        }
        return it->second;
    }

    std::pair<double, bool> getNbfix(const std::string& type1, const std::string& type2) const {
        auto key = makeTypePair(type1, type2);
        auto it = nbfix.find(key);
        if (it == nbfix.end()) {
            key = makeTypePair(type2, type1);  // Try reverse order
            it = nbfix.find(key);
            if (it == nbfix.end()) {
                return std::make_pair(0.0, false);
            }
        }
        return std::make_pair(it->second, true);
    }
};

} // namespace pygcmc


