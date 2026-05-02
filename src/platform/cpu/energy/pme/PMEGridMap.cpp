#include "PMEGridMap.hpp"
#include "PMEGlobal.hpp"
#include "PMESpline.hpp"
#include "PMECore.hpp"
#include "platform/platform.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Calculate grid indices and fractional coordinates for atoms
 */
void calculateGridIndicesAndFractions(const model::MCState& state,
                                     const std::vector<int>& atomsToProcess,
                                     const double recipBoxVectors[3][3],
                                     std::vector<std::vector<int>>& gridIndices,
                                     std::vector<std::vector<double>>& gridFractions) {
    const auto& atoms = state.atoms;
    int processedAtoms = atomsToProcess.size();

    // Resize output arrays
    gridIndices.resize(processedAtoms, std::vector<int>(3, 0));
    gridFractions.resize(processedAtoms, std::vector<double>(3, 0.0));

    // Calculate grid indices and fractional offsets for all atoms to be processed
    for (int i = 0; i < processedAtoms; i++) {
        int atomIdx = atomsToProcess[i];

        // Get position from atoms
        float pos[3] = {atoms[atomIdx].x, atoms[atomIdx].y, atoms[atomIdx].z};

        // Convert position to fractional coordinates
        double fractional[3];
        for (int d = 0; d < 3; d++) {
            // Calculate fractional coordinate - consistent with OpenMM
            // Use reciprocal lattice vectors to calculate fractional coordinates
            fractional[d] = 0.0;
            for (int j = 0; j < 3; j++) {
                fractional[d] += pos[j] * recipBoxVectors[j][d];  // No division by 2π needed
            }

            // Ensure in [0,1) range, handle periodic boundary conditions
            fractional[d] -= floor(fractional[d]);
            // Scale fractional coordinates to grid
            fractional[d] *= getPMEParams().meshSize[d];
        }

        // Calculate grid indices and fractional parts - fix: remove incorrect offset
        for (int d = 0; d < 3; d++) {
            gridFractions[i][d] = fractional[d] - floor(fractional[d]);
            // Fix: remove incorrect -order/2 offset, consistent with pme.cpp
            gridIndices[i][d] = static_cast<int>(floor(fractional[d]));
            // Ensure grid indices are within correct range
            if (gridIndices[i][d] < 0)
                gridIndices[i][d] += getPMEParams().meshSize[d];
        }
    }

    // Print grid index information for every 100 atoms - if there are enough atoms
    if (platform::is_debug_mode()) {
        for (int i = 0; i < std::min(500, processedAtoms); i += 100) {
            if (i < static_cast<int>(gridIndices.size())) {
                platform::log(LogLevel::DEBUG, "Atom " + std::to_string(atomsToProcess[i]) + " grid index: ["
                        + std::to_string(gridIndices[i][0]) + ", "
                        + std::to_string(gridIndices[i][1]) + ", "
                        + std::to_string(gridIndices[i][2]) + "]");
            }
        }
    }
}

/**
 * @brief Calculate B-spline coefficients for all atoms
 */
void calculateBSplineCoefficients(const std::vector<int>& atomsToProcess,
                                 const std::vector<std::vector<double>>& gridFractions,
                                 std::vector<std::vector<double>>& bsplines_theta) {
    int processedAtoms = atomsToProcess.size();
    int order = getPMEParams().splineOrder;

    // Initialize B-spline arrays
    bsplines_theta.resize(3);
    for (int d = 0; d < 3; d++) {
        bsplines_theta[d].resize(order * processedAtoms, 0.0);
    }

    // Calculate B-spline coefficients for all atoms to be processed
    for (int i = 0; i < processedAtoms; i++) {
        double* thetax = &bsplines_theta[0][i * order];
        double* thetay = &bsplines_theta[1][i * order];
        double* thetaz = &bsplines_theta[2][i * order];

        // Calculate B-spline coefficients for each dimension
        std::vector<double> coefficients(order);

        // X dimension B-spline
        computeBSplineCoefficients(gridFractions[i][0], order, coefficients);
        for (int j = 0; j < order; j++) {
            thetax[j] = coefficients[j];
        }

        // Y dimension B-spline
        computeBSplineCoefficients(gridFractions[i][1], order, coefficients);
        for (int j = 0; j < order; j++) {
            thetay[j] = coefficients[j];
        }

        // Z dimension B-spline
        computeBSplineCoefficients(gridFractions[i][2], order, coefficients);
        for (int j = 0; j < order; j++) {
            thetaz[j] = coefficients[j];
        }
    }
}

/**
 * @brief Distribute charges to the grid
 */
double distributeChargesToGrid(const model::MCState& state,
                              const std::vector<int>& atomsToProcess,
                              const std::vector<std::vector<int>>& gridIndices,
                              const std::vector<std::vector<double>>& bsplines_theta) {
    const auto& atoms = state.atoms;
    int processedAtoms = atomsToProcess.size();

    // Get grid dimensions
    int nx = getPMEParams().meshSize[0];
    int ny = getPMEParams().meshSize[1];
    int nz = getPMEParams().meshSize[2];
    int order = getPMEParams().splineOrder;

    // Distribute charges to grid
    double totalGridCharge = 0.0;

    // Statistics variables, only used in debug mode
    int nonZeroPoints = 0;
    int updatedPoints = 0;

    if (!platform::is_debug_mode()) {
        nonZeroPoints = -1;  // Marked as non-debug mode
        updatedPoints = -1;
    }

    for (int i = 0; i < processedAtoms; i++) {
        int atomIdx = atomsToProcess[i];
        double charge = atoms[atomIdx].charge;

        // Get grid indices and B-spline coefficients
        int x0index = gridIndices[i][0];
        int y0index = gridIndices[i][1];
        int z0index = gridIndices[i][2];

        const double* thetax = &bsplines_theta[0][i * order];
        const double* thetay = &bsplines_theta[1][i * order];
        const double* thetaz = &bsplines_theta[2][i * order];

        // Distribute charge to grid - exactly following pme_grid_spread_charge
        for (int ix = 0; ix < order; ix++) {
            int xindex = (x0index + ix) % nx;

            for (int iy = 0; iy < order; iy++) {
                int yindex = (y0index + iy) % ny;

                for (int iz = 0; iz < order; iz++) {
                    int zindex = (z0index + iz) % nz;

                    // Calculate grid index - ensure exact match with pme.cpp
                    // Original is xindex * ny * nz + yindex * nz + zindex, which is correct
                    int index = xindex * ny * nz + yindex * nz + zindex;

                    // Ensure index doesn't go out of bounds
                    if (index >= 0 && static_cast<size_t>(index) < getPMEParams().pmeGrid.size()) {
                        // Calculate B-spline weight (product of three directions)
                        double weight = thetax[ix] * thetay[iy] * thetaz[iz];

                        // Distribute charge to grid point - use same method as pme.cpp
                        double chargeContribution = charge * weight;

                        // Key fix: add contribution directly to grid point, consistent with pme.cpp
                        // In pme.cpp: pme->grid[index] += chargeContribution;
                        // This only affects the real part, as chargeContribution is real
                        getPMEParams().pmeGrid[index] += chargeContribution;

                        // Update statistics - only execute in debug mode
                        if (platform::is_debug_mode()) {
                            totalGridCharge += chargeContribution;
                            if (std::abs(chargeContribution) > 1e-10) {
                                nonZeroPoints++;

                                // Track first 10 updated grid points - modified to use same format as pme.cpp
                                if (updatedPoints < 10) {
                                    platform::log(LogLevel::DEBUG, "Updated grid point[" + std::to_string(index) + "]: charge=" + std::to_string(charge)
                                              + ", weight=" + std::to_string(weight)
                                              + ", contribution=" + std::to_string(chargeContribution));
                                    updatedPoints++;
                                }
                            }
                        } else {
                            totalGridCharge += chargeContribution; // Total charge still needs to be calculated
                        }
                    }
                }
            }
        }
    }

    // Output processing progress
    int atomIdx = state.activeAtomCount - 1; // Use index of last processed atom
    if (platform::is_debug_mode() && (atomIdx + 1) % 1000000 == 0) {
        platform::log(LogLevel::DEBUG, "Processed ", atomIdx + 1, " atoms");
    }

    // Output charge distribution completion info
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::INFO, "Charge spreading complete: ", state.activeAtomCount,
                    " atoms processed, total grid charge = ", totalGridCharge,
                    ", non-zero grid points = ", nonZeroPoints);

        platform::log(LogLevel::DEBUG, "Charge spreading complete: " + std::to_string(state.activeAtomCount)
                    + " atoms processed, total grid charge = " + std::to_string(totalGridCharge)
                    + ", non-zero grid points = " + std::to_string(nonZeroPoints));
    } else {
        // In non-debug mode, only output basic information
        platform::log(LogLevel::INFO, "Charge spreading complete: ", state.activeAtomCount,
                    " atoms processed, total grid charge = ", totalGridCharge);
    }

    return totalGridCharge;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc
