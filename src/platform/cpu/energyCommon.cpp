// src/platform/cpu/energyCommon.cpp

#include "energyCommon.hpp"
#include "energyDirect.hpp"
#include "energyEwald.hpp"
#include <stdexcept>

namespace pygcmc {
namespace platform {
namespace cpu {

// Constants for energy calculations
const float COULOMB = 138.935456f;
const float MIN_SAFE_DISTANCE = 0.01f;  // nm (1% of typical sigma)
const float MAX_SAFE_ENERGY = 1e6f;     // kJ/mol

// Debug flag
bool energy_debug_output = false;

/**
 * @brief Unified system energy calculation interface
 * 
 * @param state MC state
 * @param method Energy calculation method (DIRECT or EWALD)
 * @param use_cutoff Whether to use cutoff
 * @param use_pbc Whether to use periodic boundary conditions
 */
void computeSystemEnergy(model::MCState& state, 
                         EnergyMethod method,
                         bool use_cutoff, 
                         bool use_pbc) {
    // Validate necessary parameters for periodic boundary conditions
    if (use_pbc) {
        validateBox(state.info.box, use_cutoff ? state.info.cutoff : 0.0f);
    }
    
    // Choose different implementations based on calculation method
    switch (method) {
        case EnergyMethod::DIRECT:
            // Use direct calculation method
            computeSystemEnergyDirect(state, use_cutoff, use_pbc);
            break;
            
        case EnergyMethod::EWALD:
            // Use Ewald method (requires periodic boundary conditions)
            if (!use_pbc) {
                throw std::runtime_error("Ewald method requires periodic boundary conditions");
            }
            computeSystemEnergyEwald(state);
            break;
    }
}

/**
 * @brief Unified movement residue energy calculation interface
 * 
 * @param state MC state
 * @param method Energy calculation method (DIRECT or EWALD)
 * @param use_cutoff Whether to use cutoff
 * @param use_pbc Whether to use periodic boundary conditions
 */
void computeMovementEnergy(model::MCState& state, 
                          EnergyMethod method,
                          bool use_cutoff, 
                          bool use_pbc) {
    // Validate necessary parameters for periodic boundary conditions
    if (use_pbc) {
        validateBox(state.info.box, use_cutoff ? state.info.cutoff : 0.0f);
    }
    
    // Choose different implementations based on calculation method
    switch (method) {
        case EnergyMethod::DIRECT:
            // Use direct calculation method
            computeMovementEnergyDirect(state, use_cutoff, use_pbc);
            break;
            
        case EnergyMethod::EWALD:
            // Use Ewald method (requires periodic boundary conditions)
            if (!use_pbc) {
                throw std::runtime_error("Ewald method requires periodic boundary conditions");
            }
            computeMovementEnergyEwald(state);
            break;
    }
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 