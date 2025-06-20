#pragma once

#include "DirectCore.hpp"
#include "DirectSystemEnergy.hpp"
#include "DirectResidueEnergy.hpp"
#include "DirectPBCCalculation.hpp"
#include "model/montecarlo.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace direct {

/**
 * @brief DirectComposite - Main facade for direct energy calculation methods
 * 
 * This class provides a unified interface to all direct calculation methods,
 * following the same pattern as EwaldComposite and PMEComposite.
 * It acts as the main entry point for direct energy calculations.
 */
class DirectComposite {
public:
    /**
     * @brief Calculate system energy using direct method
     * 
     * @param state System state
     * @param use_cutoff Whether to use distance cutoff
     * @param use_pbc Whether to use periodic boundary conditions
     */
    static void calculateSystemEnergy(model::MCState& state, bool use_cutoff = false, bool use_pbc = false);
    
    /**
     * @brief Calculate movement residue energy using direct method
     * 
     * @param state System state
     * @param use_cutoff Whether to use distance cutoff
     * @param use_pbc Whether to use periodic boundary conditions
     */
    static void calculateMovementEnergy(model::MCState& state, bool use_cutoff = false, bool use_pbc = false);
    
    /**
     * @brief Calculate VDW-only energy using direct method
     * 
     * @param state System state
     * @param use_cutoff Whether to use distance cutoff
     * @param use_pbc Whether to use periodic boundary conditions
     */
    static void calculateVdwEnergy(model::MCState& state, bool use_cutoff = false, bool use_pbc = false);
};

} // namespace direct
} // namespace cpu
} // namespace platform
} // namespace pygcmc 