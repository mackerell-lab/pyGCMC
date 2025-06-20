#include "PMEGridCharge.hpp"
#include "PMESpline.hpp"
#include "PMECore.hpp"
#include "platform/platform.hpp"
#include <algorithm>
#include <cmath>
#include <vector>
#include <set>

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
    
    // Only execute the following code in debug_mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Spreading charges onto PME grid", (fixed_only ? " (fixed only)" : ""));
        
        // Calculate total system charge
        double totalCharge = 0.0;
        for (int i = 0; i < state.activeAtomCount; i++) {
            totalCharge += state.atoms[i].charge;
        }
        platform::log(LogLevel::DEBUG, "Total system charge: " + std::to_string(totalCharge));
        
        // Output initial grid values
        platform::log(LogLevel::DEBUG, "Initial values of the first 10 grid points:");
        for (int i = 0; i < 10 && i < static_cast<int>(pme_params.pmeGrid.size()); i++) {
            platform::log(LogLevel::DEBUG, "  Grid point[" + std::to_string(i) + "] = " + std::to_string(pme_params.pmeGrid[i].real()));
        }
        
        // Output charge values
        platform::log(LogLevel::DEBUG, "Charge values of the first 10 atoms:");
        for (int i = 0; i < 10 && i < state.activeAtomCount; i++) {
            platform::log(LogLevel::DEBUG, "  Atom[" + std::to_string(i) + "] charge = " + std::to_string(state.atoms[i].charge));
        }
    }
    
    // Directly access member variables instead of using getter methods
    const auto& atoms = state.atoms;
    const auto& residues = state.residues;
    // For positions, use atoms to directly access coordinates
    const auto& info = state.info;
    // Change type from double to float to match info.box type
    const float* box = info.box;
    
    // Calculate total system charge
    double totalCharge = 0.0;
    int processedAtoms = 0;
    
    // Create a list of atoms to process (either all atoms or only fixed atoms)
    std::vector<int> atomsToProcess;
    atomsToProcess.reserve(state.activeAtomCount); // Reserve space for efficiency
    
    if (fixed_only) {
        // Only include atoms from fixed residues
        for (int i = 0; i < state.activeResidueCount; ++i) {
            const auto& residue = residues[i];
            if (residue.active && residue.fixed) {
                // Add all atoms in this fixed residue
                for (int j = 0; j < residue.atomCount; ++j) {
                    int atomIndex = residue.atomStart + j;
                    atomsToProcess.push_back(atomIndex);
                    totalCharge += atoms[atomIndex].charge;
                }
            }
        }
        
        processedAtoms = atomsToProcess.size();
        platform::log(LogLevel::INFO, "Processing ", processedAtoms, " atoms from fixed residues, total charge: ", totalCharge);
    } else {
        // Process all atoms
        for (int i = 0; i < state.activeAtomCount; ++i) {
            atomsToProcess.push_back(i);
            totalCharge += atoms[i].charge;
        }
        
        processedAtoms = state.activeAtomCount;
        platform::log(LogLevel::INFO, "Total system charge: ", totalCharge);
        platform::log(LogLevel::DEBUG, "Total system charge: " + std::to_string(totalCharge));
    }
    
    // Check the initial values of the first 10 grid points
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Initial values of the first 10 grid points:");
        for (int i = 0; i < 10 && i < static_cast<int>(pme_params.pmeGrid.size()); i++) {
            platform::log(LogLevel::DEBUG, "  Grid point[" + std::to_string(i) + "] = " + std::to_string(pme_params.pmeGrid[i].real()));
        }
        
        // Print some atom charge values to verify if there are non-zero charges
        int numToPrint = std::min(10, processedAtoms);
        platform::log(LogLevel::DEBUG, "Charge values of the first " + std::to_string(numToPrint) + " processed atoms:");
        for (int i = 0; i < numToPrint; i++) {
            int atomIndex = atomsToProcess[i];
            platform::log(LogLevel::DEBUG, "  Atom[" + std::to_string(atomIndex) + "] charge = " + std::to_string(atoms[atomIndex].charge));
        }
    }
    
    // Calculate reciprocal lattice vectors - ensure consistent with pme.cpp
    double recipBoxVectors[3][3] = {{0}};
    
    // Handle periodic box vectors - simplify for diagonal boxes
    double periodicBoxVectors[3][3] = {
        {box[0], 0.0, 0.0},
        {0.0, box[1], 0.0},
        {0.0, 0.0, box[2]}
    };
    
    // Check if it's a diagonal box - no need to output
    bool isDiagonalBox = true;  // Always true since we force it to be diagonal
    
    if (isDiagonalBox) {
        // Simple calculations for diagonal boxes
        recipBoxVectors[0][0] = 2.0 * M_PI / box[0]; // 2π/a
        recipBoxVectors[1][1] = 2.0 * M_PI / box[1]; // 2π/b 
        recipBoxVectors[2][2] = 2.0 * M_PI / box[2]; // 2π/c
    } else {
        // Non-diagonal boxes require full calculation of reciprocal lattice vectors
        double det = periodicBoxVectors[0][0] * (periodicBoxVectors[1][1] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][1]) -
                     periodicBoxVectors[0][1] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][0]) +
                     periodicBoxVectors[0][2] * (periodicBoxVectors[1][0] * periodicBoxVectors[2][1] - periodicBoxVectors[1][1] * periodicBoxVectors[2][0]);
        
        // Calculate cross products and reciprocal lattice vectors
        recipBoxVectors[0][0] = 2.0 * M_PI * (periodicBoxVectors[1][1] * periodicBoxVectors[2][2] - periodicBoxVectors[1][2] * periodicBoxVectors[2][1]) / det;
        recipBoxVectors[0][1] = 2.0 * M_PI * (periodicBoxVectors[0][2] * periodicBoxVectors[2][1] - periodicBoxVectors[0][1] * periodicBoxVectors[2][2]) / det;
        recipBoxVectors[0][2] = 2.0 * M_PI * (periodicBoxVectors[0][1] * periodicBoxVectors[1][2] - periodicBoxVectors[0][2] * periodicBoxVectors[1][1]) / det;
        
        recipBoxVectors[1][0] = 2.0 * M_PI * (periodicBoxVectors[1][2] * periodicBoxVectors[2][0] - periodicBoxVectors[1][0] * periodicBoxVectors[2][2]) / det;
        recipBoxVectors[1][1] = 2.0 * M_PI * (periodicBoxVectors[0][0] * periodicBoxVectors[2][2] - periodicBoxVectors[0][2] * periodicBoxVectors[2][0]) / det;
        recipBoxVectors[1][2] = 2.0 * M_PI * (periodicBoxVectors[0][2] * periodicBoxVectors[1][0] - periodicBoxVectors[0][0] * periodicBoxVectors[1][2]) / det;
        
        recipBoxVectors[2][0] = 2.0 * M_PI * (periodicBoxVectors[1][0] * periodicBoxVectors[2][1] - periodicBoxVectors[1][1] * periodicBoxVectors[2][0]) / det;
        recipBoxVectors[2][1] = 2.0 * M_PI * (periodicBoxVectors[0][1] * periodicBoxVectors[2][0] - periodicBoxVectors[0][0] * periodicBoxVectors[2][1]) / det;
        recipBoxVectors[2][2] = 2.0 * M_PI * (periodicBoxVectors[0][0] * periodicBoxVectors[1][1] - periodicBoxVectors[0][1] * periodicBoxVectors[1][0]) / det;
    }
    
    // Reset grid - ensure all points are initialized to 0
    std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    // Get grid dimensions
    int nx = pme_params.meshSize[0];
    int ny = pme_params.meshSize[1];
    int nz = pme_params.meshSize[2];
    int order = pme_params.splineOrder;
    
    // Create temporary arrays - no resizeTempSplineArrays function
    // Create arrays for B-spline coefficients
    std::vector<std::vector<double>> bsplines_theta(3);
    for (int d = 0; d < 3; d++) {
        bsplines_theta[d].resize(order * processedAtoms, 0.0);
    }
    
    // Create arrays for grid indices and fractional parts
    std::vector<std::vector<int>> gridIndices(processedAtoms, std::vector<int>(3, 0));
    std::vector<std::vector<double>> gridFractions(processedAtoms, std::vector<double>(3, 0.0));
    
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
    
    // Calculate grid indices and fractional offsets for all atoms to be processed
    for (int i = 0; i < processedAtoms; i++) {
        int atomIdx = atomsToProcess[i];
        
        // Get position from atoms
        float pos[3] = {atoms[atomIdx].x, atoms[atomIdx].y, atoms[atomIdx].z};
        
        // Convert position to fractional coordinates
        double fractional[3];
        for (int d = 0; d < 3; d++) {
            // Calculate fractional coordinate - fix: ensure consistent with pme.cpp
            // Use reciprocal lattice vectors to calculate fractional coordinates, not simple division
            fractional[d] = 0.0;
            for (int j = 0; j < 3; j++) {
                fractional[d] += pos[j] * recipBoxVectors[j][d] / (2.0 * M_PI);
            }
            
            // Ensure in [0,1) range, handle periodic boundary conditions
            fractional[d] -= floor(fractional[d]);
            // Scale fractional coordinates to grid
            fractional[d] *= pme_params.meshSize[d];
        }
        
        // Calculate grid indices and fractional parts - fix: remove incorrect offset
        for (int d = 0; d < 3; d++) {
            gridFractions[i][d] = fractional[d] - floor(fractional[d]);
            // Fix: remove incorrect -order/2 offset, consistent with pme.cpp
            gridIndices[i][d] = static_cast<int>(floor(fractional[d]));
            // Ensure grid indices are within correct range
            if (gridIndices[i][d] < 0) 
                gridIndices[i][d] += pme_params.meshSize[d];
        }
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
        
        double* thetax = &bsplines_theta[0][i * order];
        double* thetay = &bsplines_theta[1][i * order];
        double* thetaz = &bsplines_theta[2][i * order];
        
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
                    if (index >= 0 && static_cast<size_t>(index) < pme_params.pmeGrid.size()) {
                        // Calculate B-spline weight (product of three directions)
                        double weight = thetax[ix] * thetay[iy] * thetaz[iz];
                        
                        // Distribute charge to grid point - use same method as pme.cpp
                        double chargeContribution = charge * weight;
                        
                        // Key fix: add contribution directly to grid point, consistent with pme.cpp
                        // In pme.cpp: pme->grid[index] += chargeContribution;
                        // This only affects the real part, as chargeContribution is real
                        pme_params.pmeGrid[index] += chargeContribution;
                        
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