// modules/core/include/pygcmc/core/force.hpp

#ifndef PYGCMC_CORE_FORCE_HPP
#define PYGCMC_CORE_FORCE_HPP

#include <vector>
#include <array>

namespace pygcmc {
namespace core {

// Forward declaration
class System;

/**
 * @brief Base class for all forces in the system
 */
class Force {
public:
    virtual ~Force() = default;
    
    /**
     * @brief Calculate the force and energy for the system
     * @param system The molecular system
     * @param forces Output forces array [residue_idx][atom_idx][xyz]
     * @return The potential energy of the system
     */
    virtual double calculate_forces(const System& system,
                                  std::vector<std::vector<std::array<double, 3>>>& forces) const = 0;

    /**
     * @brief Whether this force uses periodic boundary conditions
     */
    virtual bool uses_periodic_boundary_conditions() const = 0;
};

} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_FORCE_HPP 