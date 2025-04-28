// src/platform/cpu/energyPGP.hpp

#pragma once

#include "model/montecarlo.hpp"
#include "platform/platform.hpp"
#include "energyCommon.hpp"
#include "energyPME.hpp"
#include <array>
#include <vector>
#include <complex>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief PGP-PME algorithm parameter structure (Precomputed Grid-Potential Particle Mesh Ewald)
 * 
 * This structure contains all parameters and data structures required by the Precomputed Grid-Potential Particle Mesh Ewald algorithm.
 * PGP-PME is an optimized PME method that accelerates energy evaluation in Monte Carlo simulations by precomputing potential grids.
 * 
 * Main parameter categories include:
 * 1. Standard PME parameters: alpha, grid size, interpolation order, etc.
 * 2. Precomputed grid parameters: potential cutoff distance, grid size and spacing, etc.
 * 3. Data storage: precomputed potential grid, B-spline moduli, etc.
 */
struct PGPParams : public PMEParams {
    bool initialized = false;      // Whether parameters have been initialized
    double alpha;                  // Ewald separation parameter, balances real space and reciprocal space calculations
    double tolerance;              // Error tolerance
    double cutoff;                 // Real space cutoff distance
    double epsilon_r;              // Relative dielectric constant
    int splineOrder;               // B-spline interpolation order (typically 4, cubic B-spline)
    std::array<double, 3> box;     // Simulation box dimensions
    std::array<int, 3> meshSize;   // PME grid dimensions
    
    // Precomputed potential grid parameters
    double potential_cutoff;              // Cutoff distance for potential calculation
    int potential_grid_size[3];           // Precomputed potential grid dimensions
    double grid_spacing;                  // Grid spacing
    std::vector<std::complex<double>> potentialGrid;  // Precomputed potential grid data
    
    // PME algorithm parameters (mainly set by setPMEParameters function)
    std::vector<double> erfcTable;         // erfc function lookup table
    std::vector<double> ewaldScaleTable;   // Ewald scaling factor lookup table
    double ewaldDX;                        // Ewald table step size
    double ewaldDXInv;                     // Inverse of Ewald table step size
    double erfcDXInv;                      // Inverse of erfc table step size
    std::vector<double> bsplineModuli[3];  // B-spline moduli
    std::vector<std::complex<double>> pmeGrid;   // PME grid
    std::vector<double> pmeCharge;        // PME charge grid, type must match PME struct
    
    // Debug flags
    bool debug_mode = true;  // Debug mode enabled by default

    /**
     * @brief Initialize the 3D grid for precomputed potential
     * 
     * This method creates and initializes a 3D grid structure for storing precomputed potentials. The grid size is determined by the 
     * potential_grid_size parameter, typically set based on required accuracy and computational resources. Grid spacing is automatically 
     * calculated to ensure sufficient resolution in all dimensions.
     * 
     * Algorithm principles:
     * 1. Calculate total grid points based on specified grid dimensions
     * 2. Allocate memory space for storing precomputed potential values
     * 3. Calculate grid spacing, ensuring at least the required resolution in all dimensions
     * 
     * Computational complexity:
     * - Space complexity: O(nx*ny*nz), where nx/ny/nz are the grid sizes in each dimension
     * - Time complexity: O(1)
     * 
     * Usage scenarios:
     * - Automatically called after setting PGP parameters
     * - Need to reinitialize when system size changes
     */
    void initializePotentialGrid();
};

// Global PGP parameters
extern PGPParams pgp_params;

/**
 * @brief Set all parameters for the PGP-PME (Precomputed Grid-Potential Particle Mesh Ewald) algorithm
 * 
 * This function is the entry point for PGP-PME algorithm configuration, setting all parameters required to run PGP-PME.
 * It first configures standard PME parameters, then adds PGP-specific grid and cutoff parameters, and finally initializes the precomputed grid.
 * 
 * Algorithm principles:
 * 1. Call setPMEParameters to set basic PME parameters
 * 2. Copy PME parameters to the PGP parameter structure
 * 3. Add PGP-specific parameters such as potential cutoff and grid size
 * 4. Initialize the precomputed potential grid structure
 * 
 * Parameters:
 * @param alpha Ewald separation parameter, controls the balance between real space and reciprocal space calculations, typical value 0.2-0.3 Å^-1
 * @param meshSize PME calculation grid size, array form [nx,ny,nz], usually proportional to box size
 * @param potential_cutoff Cutoff distance for potential calculation, typically less than or equal to PME's real space cutoff
 * @param potentialGridSize Precomputed potential grid size, array form [nx,ny,nz], determines interpolation accuracy
 * @param splineOrder B-spline interpolation order, typically 4 (cubic B-spline), affects accuracy and calculation speed
 * @param tolerance Calculation accuracy tolerance, used for parameter selection optimization
 * 
 * Usage scenarios:
 * - Call this function for initial setup before starting simulation
 * - Reconfigure when simulation conditions (e.g., box size, accuracy requirements) change
 * - Set before starting a new Monte Carlo simulation
 */
void setPGPParameters(double alpha, const int meshSize[3], double potential_cutoff, 
                        const int potentialGridSize[3], int splineOrder, double tolerance);

/**
 * @brief Precompute the grid potential for the fixed part of the system
 * 
 * This is one of the core functions of the Precomputed Grid-Potential Particle Mesh Ewald algorithm,
 * responsible for precomputing the electrostatic potential field of the fixed part of the system. It distributes the
 * charges of the fixed part onto a grid, calculates the potential through FFT transformation, and stores the results
 * for subsequent energy calculations.
 * The precomputation step only needs to be performed once when the fixed part of the system changes, greatly
 * improving the efficiency of Monte Carlo simulations.
 * 
 * Algorithm principles:
 * 1. Distribute fixed part point charges onto a grid using B-spline interpolation
 * 2. Perform forward FFT transformation of the charge grid into reciprocal space
 * 3. Apply Ewald factors in reciprocal space for long-range corrections
 * 4. Perform reverse FFT to obtain potential distribution in real space
 * 5. Store the results in the precomputed potential grid
 * 
 * Computational complexity:
 * - Charge assignment: O(N*p^3), where N is the number of atoms, p is the B-spline order
 * - FFT: O(M*log(M)), where M is the total number of grid points (nx*ny*nz)
 * - Applying Ewald factors: O(M)
 * 
 * @param state System state, containing atom coordinates, charges, and box information
 * @param fixed_only Whether to process only the fixed part of the system (true), or all parts (false)
 * 
 * Usage scenarios:
 * - System initialization for precomputing fixed part potential
 * - In GCMC simulation, potential field of fixed part (e.g., protein) can be precomputed
 * - Need to re-call this function to update potential field when fixed part configuration changes
 */
void precomputeGridPotential(model::MCState& state, bool fixed_only = true);

/**
 * @brief Calculate the energy of moving molecules by interpolation
 * 
 * This function is another core function of the Precomputed Grid-Potential Particle Mesh Ewald algorithm,
 * used to rapidly evaluate the energy of moving molecules in the precomputed potential field. Using the precomputed
 * potential grid, it efficiently calculates the energy of moving molecules in this potential field through B-spline
 * interpolation methods, avoiding direct calculation of intermolecular interactions, greatly accelerating energy
 * evaluation in Monte Carlo simulations.
 * 
 * Algorithm principles:
 * 1. Iterate through all residues and atoms marked as moving
 * 2. For each charged atom, obtain the potential value at its position through B-spline interpolation from the precomputed potential grid
 * 3. Multiply the potential value by the atom charge and accumulate to get the total energy
 * 
 * Computational complexity:
 * - O(M*p^3), where M is the number of moving atoms, p is the B-spline order
 * - Significantly reduced compared to traditional method O(M*N), where N is the number of fixed atoms (typically N >> M)
 * 
 * @param state System state, containing information about moving molecules and the precomputed potential grid
 * @param energy Output parameter, stores calculated energy value
 * 
 * Usage scenarios:
 * - GCMC simulation for energy evaluation of molecule insertion/deletion
 * - CBMC simulation for energy comparison of different configurations
 * - MC move trial for energy change evaluation of new configuration
 */
void interpolateMoleculeEnergy(model::MCState& state, double& energy);

/**
 * @brief Calculate the energy of moving molecules by interpolation and return the result
 * 
 * This is a wrapper function for interpolateMoleculeEnergy that directly returns the calculated energy value,
 * facilitating Python calling and testing. This function internally creates an energy variable and calls the original interpolateMoleculeEnergy function.
 * 
 * @param state System state, containing information about moving molecules and the precomputed potential grid
 * @return Calculated energy value
 */
double calculateMoleculeEnergy(model::MCState& state);

/**
 * @brief Calculate real-space part of the PGP method for short-range electrostatics
 * 
 * This function calculates the short-range electrostatic interactions between atom pairs
 * using the erfc(αr)/r term, similar to the PME method's real-space component. It handles
 * the part of electrostatic interactions not covered by the precomputed grid potential.
 * 
 * @param state System state
 * @param movement_only Whether to calculate only for moving residues
 * @param store_in_residues Whether to store energy in residues
 */
void computeRealSpacePGP(model::MCState& state, bool movement_only, bool store_in_residues = true);

/**
 * @brief Calculate self energy correction for PGP method
 * 
 * This function calculates the self energy correction term in the PGP-PME method,
 * which compensates for the self-interaction that occurs in reciprocal space calculations.
 * 
 * @param state System state
 * @param movement_only Whether to calculate only for moving residues
 * @return Self energy correction value
 */
double computeSelfEnergyPGP(model::MCState& state, bool movement_only);

/**
 * @brief Use PGP method to calculate system energy
 * 
 * This function provides a complete energy calculation for the entire system using the
 * PGP-PME method, including grid potential interpolation, real-space electrostatics,
 * self energy correction, and Lennard-Jones interactions.
 * 
 * @param state MC state
 */
void computeSystemEnergyPGP(model::MCState& state);

/**
 * @brief Use PGP method to calculate energy of moving residues
 * 
 * This function calculates the energy for just the moving residues using the PGP-PME method,
 * which is useful for Monte Carlo move acceptance/rejection decisions.
 * 
 * @param state MC state
 */
void computeMovementEnergyPGP(model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc 