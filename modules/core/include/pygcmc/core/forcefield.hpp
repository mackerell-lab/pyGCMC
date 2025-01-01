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

    // Print methods
    void print_nonbonded_params() const;

    // Getters and setters
    double get_cutoff() const { return cutoff_; }
    void set_cutoff(double cutoff) { cutoff_ = cutoff; }
    
    double get_switching() const { return switching_; }
    void set_switching(double switching) { switching_ = switching; }
    
    double get_pairlist_distance() const { return pairlist_distance_; }
    void set_pairlist_distance(double distance) { pairlist_distance_ = distance; }

    // Access to parameters
    const std::map<std::string, io::ForceFieldPair>& nonbonded_params() const { return nonbonded_params_; }
    std::map<std::string, io::ForceFieldPair>& nonbonded_params() { return nonbonded_params_; }
    
    const std::map<std::pair<std::string, std::string>, io::ForceFieldPair>& nbfix_params() const { return nbfix_params_; }
    std::map<std::pair<std::string, std::string>, io::ForceFieldPair>& nbfix_params() { return nbfix_params_; }

private:
    double cutoff_ = 14.0;
    double switching_ = 12.0;
    double pairlist_distance_ = 16.0;
    std::map<std::string, io::ForceFieldPair> nonbonded_params_;
    std::map<std::pair<std::string, std::string>, io::ForceFieldPair> nbfix_params_;
}; 

} // namespace core
} // namespace pygcmc 