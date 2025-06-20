#include "EnergySystemInterface.hpp"
#include "../direct/DirectComposite.hpp"
#include "../direct/DirectPBCCalculation.hpp"
#include "../direct/DirectSystemEnergy.hpp"
#include "../ewald/EwaldComposite.hpp"  // Include inline function definitions
#include "../pme/PMEComposite.hpp"      // Include inline function definitions
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
            // Use direct calculation method via DirectComposite
            direct::DirectComposite::calculateSystemEnergy(state, use_cutoff, use_pbc);
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
            // Use direct calculation method via DirectComposite
            direct::DirectComposite::calculateMovementEnergy(state, use_cutoff, use_pbc);
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

/**
 * @brief Forward to direct::computeSystemEnergyPBC
 */
void computeSystemEnergyPBC(model::MCState& state) {
    direct::computeSystemEnergyPBC(state);
}

/**
 * @brief Forward to direct::computeSystemEnergyPBCCutoff
 */
void computeSystemEnergyPBCCutoff(model::MCState& state) {
    direct::computeSystemEnergyPBCCutoff(state);
}

/**
 * @brief Forward to direct::computeSystemVdwEnergyCutoff
 */
void computeSystemVdwEnergyCutoff(model::MCState& state) {
    direct::computeSystemVdwEnergyCutoff(state);
}

/**
 * @brief Forward to direct::computeSystemVdwEnergyDirect
 */
void computeSystemVdwEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc) {
    direct::computeSystemVdwEnergyDirect(state, use_cutoff, use_pbc);
}

/**
 * @brief Legacy function - use unified interface instead
 */
void computeMovementEnergyCutoff(model::MCState& state) {
    computeMovementEnergy(state, EnergyMethod::DIRECT, true, false);
}

/**
 * @brief Legacy function - use unified interface instead
 */
void computeSystemEnergyCutoff(model::MCState& state) {
    computeSystemEnergy(state, EnergyMethod::DIRECT, true, false);
}

// Note: Other energy functions are already defined in their respective modules
// This file only implements the unified interface and missing direct calculation functions

} // namespace cpu
} // namespace platform
} // namespace pygcmc 