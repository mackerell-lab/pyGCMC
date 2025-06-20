#include "EnergyInterface.hpp"
#include "EnergyDirectCore.hpp"
#include "../ewald/EwaldMain.hpp"  // Include inline function definitions
#include "../pme/PMEMain.hpp"      // Include inline function definitions
#include "EnergyUtils.hpp"
#include <stdexcept>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Unified system energy calculation interface
 * 
 * @param state MC state
 * @param method Energy calculation method (DIRECT, EWALD, or PME)
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
            
        case EnergyMethod::PME:
            // Use Particle Mesh Ewald method (requires periodic boundary conditions)
            if (!use_pbc) {
                throw std::runtime_error("PME method requires periodic boundary conditions");
            }
            computeSystemEnergyPME(state);
            break;
    }
}

/**
 * @brief Unified movement residue energy calculation interface
 * 
 * @param state MC state
 * @param method Energy calculation method (DIRECT, EWALD, or PME)
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
            
        case EnergyMethod::PME:
            // Use Particle Mesh Ewald method (requires periodic boundary conditions)
            if (!use_pbc) {
                throw std::runtime_error("PME method requires periodic boundary conditions");
            }
            computeMovementEnergyPME(state);
            break;
    }
}

// Note: Direct summation functions are now directly available from DirectSummation.hpp
// No forwarding needed as they are in the same namespace

// Note: Legacy functions like computeMovementEnergyCutoff and computeSystemEnergyCutoff
// are now implemented in DirectSummation.cpp to avoid duplication

} // namespace cpu
} // namespace platform
} // namespace pygcmc 