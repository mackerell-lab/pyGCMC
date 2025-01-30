// src/model/forcefield.hpp

#pragma once

#include <string>
#include <map>
#include <vector>
#include <tuple>
#include <stdexcept>

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
};

} // namespace pygcmc


