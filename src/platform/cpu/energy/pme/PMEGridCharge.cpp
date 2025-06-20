#include "PMEGridCharge.hpp"
#include "PMEGridPrepare.hpp"
#include "PMEGridMapping.hpp"
#include "PMECore.hpp"
#include "platform/platform.hpp"
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Spread charges onto the PME grid
 * 
 * @param state MC state
 * @param fixed_only If true, only consider charges from atoms in fixed residues
 */
void spreadChargesOntoGrid(model::MCState& state, bool fixed_only) {
    // Log the start of processing
    platform::log(LogLevel::DEBUG, "Spreading charges onto PME grid", (fixed_only ? " (fixed only)" : ""));
    
    // Reset grid - ensure all points are initialized to 0
    std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    // Step 1: Select atoms to process and calculate total charge
    std::vector<int> atomsToProcess;
    double totalCharge = 0.0;
    selectAtomsToProcess(state, fixed_only, atomsToProcess, totalCharge);
    
    // Step 2: Output debug information
    outputDebugInfo(state, atomsToProcess);
    
    // Step 3: Calculate reciprocal lattice vectors
    const float* box = state.info.box;
    double recipBoxVectors[3][3];
    calculateReciprocalLatticeVectors(box, recipBoxVectors);
    
    // Step 4: Calculate grid indices and fractional coordinates
    std::vector<std::vector<int>> gridIndices;
    std::vector<std::vector<double>> gridFractions;
    calculateGridIndicesAndFractions(state, atomsToProcess, recipBoxVectors, gridIndices, gridFractions);
    
    // Step 5: Calculate B-spline coefficients
    std::vector<std::vector<double>> bsplines_theta;
    calculateBSplineCoefficients(atomsToProcess, gridFractions, bsplines_theta);
    
    // Step 6: Distribute charges to grid
    double totalGridCharge = distributeChargesToGrid(state, atomsToProcess, gridIndices, bsplines_theta);
    
    // Suppress unused variable warning
    (void)totalGridCharge;
}

/**
 * @brief Calculate fractional coordinates in the grid
 */
void calculateFractionalCoordinates(const float position[3], 
                                  const float box[3], 
                                  const int meshSize[3], 
                                  double fractional[3]) {
    for (int d = 0; d < 3; d++) {
        fractional[d] = position[d] / box[d];
        // Ensure in [0,1) range
        fractional[d] -= floor(fractional[d]);
        // Scale to grid
        fractional[d] *= meshSize[d];
    }
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 