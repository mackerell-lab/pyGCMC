// src/model/forcefield.hpp

#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <array>
#include <stdexcept>
#include <tuple>

namespace pygcmc {

struct NonbondedParams {
    int nbxmod = 5;           // interaction modification flag
    bool cdiel = false;       // constant dielectric
    bool fshift = false;      // force shifting
    bool vatom = false;          
    bool vdistance = false;      
    bool vfswitch = false;       
    double cutnb = 14.0;     // nonbonded cutoff
    double ctofnb = 12.0;    // outer cutoff for switching
    double ctonnb = 10.0;    // inner cutoff for switching
    double eps = 1.0;        // dielectric constant
    double e14fac = 1.0;     // 1-4 interaction scaling factor
    double wmin = 1.5;       // minimum weight for switching
};

struct LJParams {
    double epsilon = 0.0;    // well depth
    double rmin = 0.0;       // minimum energy distance
};

class ForceField {
public:
    ForceField() = default;
    ~ForceField() = default;

    // Parse parameters from .str or .prm file
    void parseFromFile(const std::string& filename);

    // Getters with Python-style naming
    NonbondedParams& get_nonbonded_params() { return nonbondedParams_; }
    const NonbondedParams& get_nonbonded_params() const { return nonbondedParams_; }
    
    const LJParams& get_lj_params(const std::string& atomType) const {
        auto it = ljParams_.find(atomType);
        if (it == ljParams_.end()) {
            throw std::runtime_error("LJ parameters not found for atom type: " + atomType);
        }
        return it->second;
    }

    std::tuple<double, bool> get_nbfix(const std::string& type1, const std::string& type2) const {
        auto it1 = nbfix_.find(type1);
        if (it1 != nbfix_.end()) {
            auto it2 = it1->second.find(type2);
            if (it2 != it1->second.end()) {
                return std::make_tuple(it2->second[0], true);  // return epsilon and found
            }
        }
        return std::make_tuple(0.0, false);
    }

    // Setters with Python-style naming
    void add_lj_params(const std::string& atomType, const LJParams& params) {
        ljParams_[atomType] = params;
    }

    void add_nbfix(const std::string& type1, const std::string& type2, double epsilon, double rmin) {
        nbfix_[type1][type2] = {epsilon, rmin};
        nbfix_[type2][type1] = {epsilon, rmin};  // Add symmetric pair
    }

private:
    NonbondedParams nonbondedParams_;
    std::unordered_map<std::string, LJParams> ljParams_;        // atom type -> LJ parameters
    std::unordered_map<std::string, std::unordered_map<std::string, std::array<double, 2>>> nbfix_;  // type1->type2->[epsilon, rmin]

    // Helper functions for parsing
    void parseNonbondedSection(std::istream& input);
    void parseNBFixSection(std::istream& input);
};

} // namespace pygcmc


