#include "pygcmc/core/forcefield.hpp"
#include <iostream>
#include <iomanip>

namespace pygcmc {
namespace core {

void ForceField::print_nonbonded_params() const {
    std::cout << "\nNonbonded Parameters Loaded:" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    std::cout << std::left << std::setw(10) << "Type" 
              << std::setw(15) << "Epsilon" 
              << std::setw(15) << "Rmin" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    
    for (const auto& [atomType, params] : nonbonded_params_) {
        std::cout << std::left << std::setw(10) << atomType 
                  << std::setw(15) << params.epsilon 
                  << std::setw(15) << params.rmin << std::endl;
    }
    std::cout << std::endl;
}

} // namespace core
} // namespace pygcmc 