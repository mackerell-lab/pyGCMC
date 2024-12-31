// modules/core/include/pygcmc/core/forcefield.hpp

#pragma once

#include <string>
#include <map>
#include <memory>
#include "pygcmc/core/io/ff_parser.hpp"

namespace pygcmc {
namespace core {

class ForceField {
public:
    ForceField() = default;
    ~ForceField() = default;

    // Global parameters
    double get_cutoff() const { return cutoff_; }
    double get_switching() const { return switching_; }
    double get_pairlist_distance() const { return pairlist_distance_; }
    void set_cutoff(double value) { cutoff_ = value; }
    void set_switching(double value) { switching_ = value; }
    void set_pairlist_distance(double value) { pairlist_distance_ = value; }

    // Parameter access
    const std::map<std::string, io::ForceFieldPair>& get_nonbonded_params() const { 
        return nonbonded_params_; 
    }
    const std::map<std::pair<std::string, std::string>, io::ForceFieldPair>& get_nbfix_params() const { 
        return nbfix_params_; 
    }

    // Python interface helpers
    std::map<std::string, io::ForceFieldPair>& nonbonded_params() { return nonbonded_params_; }
    std::map<std::pair<std::string, std::string>, io::ForceFieldPair>& nbfix_params() { 
        return nbfix_params_; 
    }

private:
    double cutoff_ = 14.0;
    double switching_ = 12.0;
    double pairlist_distance_ = 16.0;
    std::map<std::string, io::ForceFieldPair> nonbonded_params_;
    std::map<std::pair<std::string, std::string>, io::ForceFieldPair> nbfix_params_;
}; 

} // namespace core
} // namespace pygcmc 