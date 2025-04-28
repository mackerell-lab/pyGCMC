// src/platform/cpu/energyPGP.cpp

#include "energyPGP.hpp"
#include "energyPME.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

/**
 * @file energyPGP.cpp
 * @brief Implementation of Precomputed Grid-Potential Particle Mesh Ewald for Monte Carlo (PGP-PME-MC)
 * 
 * PGP-PME (Precomputed Grid-Potential Particle Mesh Ewald) is a particle mesh Ewald method optimized for Monte Carlo simulations.
 * It greatly accelerates the calculation of energy changes for small molecules in MC simulations by precomputing the electrostatic potential grid of the fixed parts of the system.
 * 
 * Key features of the algorithm:
 * 1. Precomputed: The electrostatic potential of the fixed parts of the system is precomputed and stored in a grid
 * 2. Grid-Potential: Uses a three-dimensional grid to represent potential distribution, combined with B-spline interpolation for efficient sampling
 * 3. PME-based: Based on traditional PME method for handling long-range electrostatic interactions, but optimized for MC
 * 
 * Workflow:
 * - Initialization phase: Set up regular PME parameters and additional potential grid parameters
 * - Precomputation phase: Distribute fixed part charges to the grid, perform FFT and store the grid potential
 * - MC simulation phase: Quickly evaluate energy changes of moving molecules through grid interpolation
 * 
 * Application scenarios:
 * - GCMC (Grand Canonical Monte Carlo) simulations for solvent molecule insertion/deletion
 * - CBMC (Configurational Bias Monte Carlo) for conformational sampling
 * - Biomolecular systems such as protein-ligand interaction simulations
 */

namespace pygcmc {
namespace platform {
namespace cpu {

// Initialize global PGP parameters
PGPParams pgp_params;

/**
 * @brief Initialize the three-dimensional grid for precomputed potential
 * 
 * This is a fundamental step in the PGP-PME algorithm, responsible for creating and initializing the 3D grid used to store precomputed potentials.
 * This function allocates grid memory based on potential_grid_size parameters and calculates appropriate grid spacing.
 */
void PGPParams::initializePotentialGrid() {
    // Calculate total grid size and allocate memory
    int totalSize = potential_grid_size[0] * potential_grid_size[1] * potential_grid_size[2];
    potentialGrid.resize(totalSize);
    
    // Set grid spacing, taking the minimum value of the three dimensions
    grid_spacing = std::min({
        box[0] / potential_grid_size[0],
        box[1] / potential_grid_size[1],
        box[2] / potential_grid_size[2]
    });
    
    // Only output debug information in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "PGP potential grid initialized with size: ", 
                    potential_grid_size[0], "x", potential_grid_size[1], "x", potential_grid_size[2],
                    ", grid spacing: ", grid_spacing);
    }
}

/**
 * @brief Set all parameters for the PGP-PME algorithm
 * 
 * This function is the entry point for the PGP-PME algorithm, used to configure all parameters needed for the algorithm to run.
 * It first sets standard PME parameters, then adds PGP-specific parameters, and finally initializes the potential grid.
 * 
 * @param alpha Ewald separation parameter, controls the balance between real-space and reciprocal-space calculations
 * @param meshSize Regular PME grid size
 * @param potential_cutoff Cutoff distance for potential calculation
 * @param potentialGridSize Grid size for precomputed potential
 * @param splineOrder Order of B-spline interpolation
 * @param tolerance Tolerance for calculation precision
 */
void setPGPParameters(double alpha, const int meshSize[3], double potential_cutoff, 
                        const int potentialGridSize[3], int splineOrder, double tolerance) {
    // First set standard PME parameters
    setPMEParameters(alpha, meshSize, splineOrder, tolerance);
    
    // Copy standard PME parameters to PGP parameter structure
    pgp_params.alpha = pme_params.alpha;
    pgp_params.tolerance = pme_params.tolerance;
    pgp_params.initialized = pme_params.initialized;
    pgp_params.cutoff = pme_params.cutoff;
    pgp_params.epsilon_r = pme_params.epsilon_r;
    pgp_params.splineOrder = pme_params.splineOrder;
    
    // Copy box size and grid size
    for (int i = 0; i < 3; i++) {
        pgp_params.box[i] = pme_params.box[i];
        pgp_params.meshSize[i] = pme_params.meshSize[i];
    }
    
    // Copy PME lookup tables
    pgp_params.erfcTable = pme_params.erfcTable;
    pgp_params.ewaldScaleTable = pme_params.ewaldScaleTable;
    pgp_params.ewaldDX = pme_params.ewaldDX;
    pgp_params.ewaldDXInv = pme_params.ewaldDXInv;
    pgp_params.erfcDXInv = pme_params.erfcDXInv;
    
    // Copy B-spline moduli
    for (int i = 0; i < 3; i++) {
        pgp_params.bsplineModuli[i] = pme_params.bsplineModuli[i];
    }
    
    // Copy PME grid
    pgp_params.pmeGrid = pme_params.pmeGrid;
    pgp_params.pmeCharge = pme_params.pmeCharge;
    
    // Set PGP-specific parameters
    pgp_params.potential_cutoff = potential_cutoff;
    for (int i = 0; i < 3; i++) {
        pgp_params.potential_grid_size[i] = potentialGridSize[i];
    }
    
    // Initialize grid for precomputed potential
    pgp_params.initializePotentialGrid();
    
    // Only output parameter setting information in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "PGP parameters set: alpha=", alpha, 
                    ", potential_cutoff=", potential_cutoff, 
                    ", potentialGrid=[", potentialGridSize[0], ",", potentialGridSize[1], ",", potentialGridSize[2], "]");
    }
}

/**
 * @brief Precompute grid potential for fixed parts of the system
 * 
 * This is one of the core functions of the PGP-PME algorithm, responsible for precomputing the electrostatic potential field of the fixed parts of the system.
 * This function assigns charges from the fixed parts to the grid, computes the potential through FFT transformation, and stores the result for later use.
 * The precomputation step needs only to be executed when the fixed part of the system changes, significantly improving the efficiency of MC simulations.
 * 
 * @param state System state, including atom coordinates, charges, and box information
 * @param fixed_only Whether to process only the fixed parts of the system (true) or all parts (false)
 */
void precomputeGridPotential(model::MCState& state, bool fixed_only) {
    // Check if parameters are initialized
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // Only output debug information in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Starting grid potential precomputation");
        platform::log(LogLevel::DEBUG, "Processing: ", (fixed_only ? "fixed parts only" : "all parts"));
        platform::log(LogLevel::DEBUG, "Grid size: ", pgp_params.potential_grid_size[0], "x", 
                     pgp_params.potential_grid_size[1], "x", 
                     pgp_params.potential_grid_size[2]);
    }
    
    // Backup PME grid, will restore later
    std::vector<std::complex<double>> pmeGridBackup = pme_params.pmeGrid;
    
    // Reset PME grid, prepare for new calculation
    std::fill(pme_params.pmeGrid.begin(), pme_params.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    // Statistics - only calculated in debug mode
    int fixed_residues_count = 0;
    
    // Calculate number of fixed residues - only computed in debug mode or when checking fixed_only validity
    if (platform::is_debug_mode() || fixed_only) {
        for (int i = 0; i < state.activeResidueCount; ++i) {
            const auto& res = state.residues[i];
            if (res.fixed && res.active) fixed_residues_count++;
        }
        
        if (platform::is_debug_mode()) {
            platform::log(LogLevel::DEBUG, "Number of fixed residues: ", fixed_residues_count);
            
            // Print information for all residues for debugging
            platform::log(LogLevel::DEBUG, "Printing fixed status for all residues:");
            for (int i = 0; i < state.activeResidueCount; ++i) {
                const auto& res = state.residues[i];
                platform::log(LogLevel::DEBUG, "Residue ", i, ": fixed=", res.fixed, ", active=", res.active,
                            ", atomCount=", res.atomCount);
            }
        }
    }
    
    // If no fixed residues but requesting fixed_only, issue warning and switch automatically
    if (fixed_only && fixed_residues_count == 0) {
        platform::log(LogLevel::WARNING, "No fixed residues found, switching to process all atoms");
        fixed_only = false;
    }
    
    // Set pme_params grid size to match pgp_params grid size, ensuring calculation uses the same grid
    for (int i = 0; i < 3; i++) {
        pme_params.meshSize[i] = pgp_params.potential_grid_size[i];
    }
    
    // Adjust pme_params grid size to fit new grid dimensions
    int totalGridSize = pgp_params.potential_grid_size[0] * pgp_params.potential_grid_size[1] * pgp_params.potential_grid_size[2];
    pme_params.pmeGrid.resize(totalGridSize, std::complex<double>(0.0, 0.0));
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Calling PME charge spreading function (fixed_only=", fixed_only, ")");
    }
    
    spreadChargesOntoGrid(state, fixed_only);

    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Calling PME forward FFT function");
    }
    
    performFFTForward();

    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Manually applying Ewald factor");
    }
    
    // Get grid size and box size
    int nx = pgp_params.potential_grid_size[0];
    int ny = pgp_params.potential_grid_size[1];
    int nz = pgp_params.potential_grid_size[2];
    double volume = pgp_params.box[0] * pgp_params.box[1] * pgp_params.box[2];
    
    // Calculate constants needed for Ewald factor application
    double alpha = pgp_params.alpha;
    // Correction: Correct coefficient in exp(-k²/(4α²))
    double factor = 1.0/(4.0*alpha*alpha);
    
    // Only print key parameters in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Key calculation parameters:");
        platform::log(LogLevel::DEBUG, "Box size: [", pgp_params.box[0], ", ", pgp_params.box[1], ", ", pgp_params.box[2], "] nm");
        platform::log(LogLevel::DEBUG, "Box volume(Ω): ", volume, " nm³");
        platform::log(LogLevel::DEBUG, "Grid size: [", nx, ", ", ny, ", ", nz, "]");
        platform::log(LogLevel::DEBUG, "Total grid points: ", nx * ny * nz);
        platform::log(LogLevel::DEBUG, "Ewald separation parameter(α): ", alpha, " nm⁻¹");
        platform::log(LogLevel::DEBUG, "exp(-k²/(4α²)) coefficient: ", factor);
    }
    
    // Get maximum k vector index
    int maxkx = (nx+1)/2;
    int maxky = (ny+1)/2;
    int maxkz = (nz+1)/2;
    
    // Calculate reciprocal lattice vectors
    double recipBoxVectors[3][3] = {{0}};
    // Correction: Reciprocal lattice vectors should be 2π/box, not 1/box
    // PME theory defines k vectors as k = 2π·n/L, missing 2π results in m² value being smaller
    recipBoxVectors[0][0] = 2.0 * M_PI / pgp_params.box[0]; 
    recipBoxVectors[1][1] = 2.0 * M_PI / pgp_params.box[1]; 
    recipBoxVectors[2][2] = 2.0 * M_PI / pgp_params.box[2];
    
    // Only print reciprocal lattice vectors and example k point values in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Reciprocal lattice vectors: [", recipBoxVectors[0][0], ", ", recipBoxVectors[1][1], ", ", recipBoxVectors[2][2], "] nm⁻¹");
        
        // Example calculation of values for a few k points
        platform::log(LogLevel::DEBUG, "Example k point values (nx/4, ny/4, nz/4):");
        int sx = nx/4, sy = ny/4, sz = nz/4;
        double mkx = sx * recipBoxVectors[0][0];
        double mky = sy * recipBoxVectors[1][1];
        double mkz = sz * recipBoxVectors[2][2];
        double mk2 = mkx*mkx + mky*mky + mkz*mkz;
        platform::log(LogLevel::DEBUG, "k = [", mkx, ", ", mky, ", ", mkz, "] nm⁻¹");
        platform::log(LogLevel::DEBUG, "|k|² = ", mk2, " nm⁻²");
        platform::log(LogLevel::DEBUG, "exp(-k²/(4α²)) = ", exp(-mk2 * factor));
    }
    
    // Apply Ewald factor
    for (int kx = 0; kx < nx; kx++) {
        double mx = (kx < maxkx) ? kx : (kx-nx);
        double mhx = mx * recipBoxVectors[0][0];
        
        for (int ky = 0; ky < ny; ky++) {
            double my = (ky < maxky) ? ky : (ky-ny);
            double mhy = my * recipBoxVectors[1][1];
            
            for (int kz = 0; kz < nz; kz++) {
                // Skip zero frequency
                if (kx == 0 && ky == 0 && kz == 0) {
                    continue;
                }
                
                double mz = (kz < maxkz) ? kz : (kz-nz);
                double mhz = mz * recipBoxVectors[2][2];
                
                // Grid index
                int index = kx * ny * nz + ky * nz + kz;
                // Get structure factor
                std::complex<double> structureFactor = pme_params.pmeGrid[index];
                
                // Calculate |k|^2
                double m2 = mhx * mhx + mhy * mhy + mhz * mhz;
                
                // Apply B-spline coefficients
                double bx = pgp_params.bsplineModuli[0][kx];
                double by = pgp_params.bsplineModuli[1][ky];
                double bz = pgp_params.bsplineModuli[2][kz];
                double denom = m2 * bx * by * bz; // Removed boxfactor, consistent with PME theory
                
                // Avoid division by zero problem
                if (denom < 1e-10) {
                    denom = 1e-10;
                }
                
                // Only apply k-dependent factor: exp(-k²/(4α²))/(k² · B)
                // Constant factor(4π/Ω) will be applied in inverse FFT
                double kDependentFactor = exp(-m2 * factor) / denom;
                
                // Apply k-dependent factor
                pme_params.pmeGrid[index] = structureFactor * kDependentFactor;
            }
        }
    }
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Executing inverse FFT to get real space potential");
    }
    
    performFFTBackward();
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Compensating inverse FFT normalization factor and applying constant factor");
    }
    
    int totalFFTPoints = nx * ny * nz;
    // Constant factor: 4π/Ω
    double constantFactor = 4.0 * M_PI / volume;
    
    // Potential physical unit conversion factor
    double ONE_4PI_EPS0 = 138.935456; // kJ·mol^-1·nm·e^-2, consistent with PME definition
    double physicalUnitFactor = ONE_4PI_EPS0 / pgp_params.epsilon_r;
    
    // Total correction factor = FFT normalization compensation(N) × 4π/Ω constant factor × physical unit conversion × 0.5(halving)
    double totalFactor = totalFFTPoints * constantFactor * physicalUnitFactor * 0.5;
    
    // Only print correction factor information in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Applying correction factor:");
        platform::log(LogLevel::DEBUG, "FFT normalization compensation(N): ", totalFFTPoints);
        platform::log(LogLevel::DEBUG, "4π/Ω constant factor: ", constantFactor, " nm⁻³");
        platform::log(LogLevel::DEBUG, "Physical unit conversion: ", physicalUnitFactor, " kJ*nm/mol*e²");
        platform::log(LogLevel::DEBUG, "Potential halving factor: 0.5");
        platform::log(LogLevel::DEBUG, "Total correction factor: ", totalFactor);
    }
    
    // Only print pre-correction potential range in debug mode
    if (platform::is_debug_mode()) {
        double preMin = 0.0, preMax = 0.0, preSum = 0.0;
        bool firstValue = true;
        
        for (int i = 0; i < totalGridSize; i++) {
            double val = pme_params.pmeGrid[i].real();
            if (firstValue) {
                preMin = preMax = val;
                firstValue = false;
            } else {
                preMin = std::min(preMin, val);
                preMax = std::max(preMax, val);
            }
            preSum += val;
        }
        platform::log(LogLevel::DEBUG, "Pre-correction potential range: [", preMin, ", ", preMax, "], Average: ", (preSum/totalGridSize));
    }
    
    // Apply total correction factor to each grid point
    for (int i = 0; i < totalGridSize; i++) {
        pme_params.pmeGrid[i] *= totalFactor;
    }
    
    // Only print post-correction potential range in debug mode
    if (platform::is_debug_mode()) {
        double postMin = 0.0, postMax = 0.0, postSum = 0.0;
        bool firstValue = true;
        
        for (int i = 0; i < totalGridSize; i++) {
            double val = pme_params.pmeGrid[i].real();
            if (firstValue) {
                postMin = postMax = val;
                firstValue = false;
            } else {
                postMin = std::min(postMin, val);
                postMax = std::max(postMax, val);
            }
            postSum += val;
        }
        platform::log(LogLevel::DEBUG, "Post-correction potential range: [", postMin, ", ", postMax, "], Average: ", (postSum/totalGridSize));
    }
    
    // Copy modified PME grid to PGP's potentialGrid
    pgp_params.potentialGrid = pme_params.pmeGrid;

    // Restore original PME grid
    pme_params.pmeGrid = pmeGridBackup;
    
    // Verify potential grid has reasonable values - only executed in debug mode
    if (platform::is_debug_mode()) {
        int potentials_nonzero = 0;
        double max_potential = 0.0;
        double min_potential = 0.0;
        double sum_potential = 0.0;
        bool first_pot = true;
        
        // Check grid point values
        for (const auto& val : pgp_params.potentialGrid) {
            double pot_val = val.real();
            sum_potential += pot_val;
            if (std::abs(pot_val) > 1e-10) {
                potentials_nonzero++;
                if (first_pot) {
                    max_potential = min_potential = pot_val;
                    first_pot = false;
                } else {
                    max_potential = std::max(max_potential, pot_val);
                    min_potential = std::min(min_potential, pot_val);
                }
            }
        }
        
        // Output potential grid statistics
        platform::log(LogLevel::DEBUG, "Potential grid statistics:");
        platform::log(LogLevel::DEBUG, "   Non-zero points: ", potentials_nonzero);
        platform::log(LogLevel::DEBUG, "   Maximum potential: ", max_potential);
        platform::log(LogLevel::DEBUG, "   Minimum potential: ", min_potential);
        platform::log(LogLevel::DEBUG, "   Potential sum: ", sum_potential);
        
        platform::log(LogLevel::DEBUG, "Potential precomputation completed");
        platform::log(LogLevel::DEBUG, "   Non-zero points: ", potentials_nonzero);
        platform::log(LogLevel::DEBUG, "   Potential range: [", min_potential, ", ", max_potential, "]");
    }
}

/**
 * @brief Calculate moving molecule energy through interpolation
 * 
 * This function is another core function of the PGP-PME algorithm, used to quickly evaluate the energy of moving molecules in the precomputed potential field.
 * By using B-spline interpolation from the precomputed grid potential, it avoids direct calculation of intermolecular interactions,
 * greatly improving the efficiency of energy evaluation in MC simulations.
 * 
 * @param state System state, containing information about moving molecules and precomputed potential grid
 * @param energy Output parameter that stores the calculated energy value
 */
void interpolateMoleculeEnergy(model::MCState& state, double& energy) {
    // Check if parameters are initialized
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized");
    }
    
    // Only output logs in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Calculate moving molecule energy through interpolation");
        platform::log(LogLevel::DEBUG, "Number of movement residue groups: ", state.movementResidues.size());
    }
    
    // Reset energy accumulator
    energy = 0.0;
    
    // Check if precomputed grid is empty - only perform full check in debug mode
    bool gridEmpty = false;
    if (platform::is_debug_mode()) {
        gridEmpty = true;
        for (const auto& val : pgp_params.potentialGrid) {
            if (std::abs(val.real()) > 1e-10 || std::abs(val.imag()) > 1e-10) {
                gridEmpty = false;
                break;
            }
        }
        
        if (gridEmpty) {
            platform::log(LogLevel::WARNING, "PGP grid is empty or not correctly initialized!");
            
            platform::log(LogLevel::DEBUG, "Grid sample point values:");
            for (int i = 0; i < std::min(10, static_cast<int>(pgp_params.potentialGrid.size())); i++) {
                platform::log(LogLevel::DEBUG, "Grid point ", i, ": ", pgp_params.potentialGrid[i].real());
            }
        } else {
            platform::log(LogLevel::DEBUG, "PGP grid contains non-zero values");
        }
    } else {
        // Non-debug mode performs simple check
        if (!pgp_params.potentialGrid.empty() && 
            std::abs(pgp_params.potentialGrid[0].real()) < 1e-10 && 
            std::abs(pgp_params.potentialGrid[0].imag()) < 1e-10) {
            // Only check the first element as a quick check
            platform::log(LogLevel::WARNING, "PGP grid may be empty or not correctly initialized!");
        }
    }
    
    // Check movement residue settings
    if (state.movementResidues.empty()) {
        platform::log(LogLevel::WARNING, "No movement residue information set!");
        
        if (platform::is_debug_mode()) {
            platform::log(LogLevel::DEBUG, "Attempting to use non-fixed residues as movement residues");
            
            // Print information about each residue, help with debugging
            platform::log(LogLevel::DEBUG, "Residue status:");
            for (int i = 0; i < state.activeResidueCount; ++i) {
                const auto& res = state.residues[i];
                platform::log(LogLevel::DEBUG, "Residue ", i, ": fixed=", res.fixed, 
                            ", active=", res.active, 
                            ", atomCount=", res.atomCount);
            }
        }
    }
    
    // Statistics - only used for debug logs
    int totalAtoms = 0;
    int chargedAtoms = 0;
    
    double raw_energy = 0.0; // Used to store un-scaled energy
    
    // Process each movement residue in the system
    if (state.movementResidues.empty()) {
        if (platform::is_debug_mode()) {
            platform::log(LogLevel::DEBUG, "No movement residue information set!");
            platform::log(LogLevel::DEBUG, "Attempting to find all non-fixed active residues...");
        }
        
        for (size_t i = 0; i < state.residues.size(); i++) {
            const auto& residue = state.residues[i];
            
            if (!residue.active || residue.fixed) continue;
            
            if (platform::is_debug_mode()) {
                platform::log(LogLevel::DEBUG, "Using non-fixed residue ", i);
            }
            
            // Process atoms...
            for (int j = 0; j < residue.atomCount; j++) {
                int atom_index = residue.atomStart + j;
                const auto& atom = state.atoms[atom_index];
                
                // Atom count only used for debug logs
                if (platform::is_debug_mode()) {
                    totalAtoms++;
                }
                
                // Only process charged atoms
                if (std::abs(atom.charge) < 1e-6) continue;
                
                // Charged atom count only used for debug logs
                if (platform::is_debug_mode()) {
                    chargedAtoms++;
                }
                
                if (platform::is_debug_mode()) {
                    platform::log(LogLevel::DEBUG, "Processing atom ", atom_index, ": position=(", 
                                atom.x, ",", atom.y, ",", atom.z, 
                                "), charge=", atom.charge);
                }
                
                // Calculate grid position and B-spline interpolation weights
                // Get atom position
                double pos[3] = {atom.x, atom.y, atom.z};
                
                // Calculate fractional coordinates - directly use box dimensions for transformation
                double fractional[3];
                for (int d = 0; d < 3; d++) {
                    fractional[d] = pos[d] / pgp_params.box[d];
                    fractional[d] -= floor(fractional[d]);  // Ensure in [0,1) range
                    fractional[d] *= pgp_params.potential_grid_size[d]; // Scale to grid
                }
                
                // Calculate grid index and fractional part
                int gridIndices[3];
                double gridFractions[3];
                for (int d = 0; d < 3; d++) {
                    gridFractions[d] = fractional[d] - floor(fractional[d]);
                    gridIndices[d] = static_cast<int>(floor(fractional[d]));
                    // Ensure grid index within correct range
                    if (gridIndices[d] < 0) 
                        gridIndices[d] += pgp_params.potential_grid_size[d];
                }
                
                // Calculate B-spline coefficients
                int nx = pgp_params.potential_grid_size[0];
                int ny = pgp_params.potential_grid_size[1];
                int nz = pgp_params.potential_grid_size[2];
                int order = pgp_params.splineOrder;
                
                std::vector<double> thetaX(order);
                std::vector<double> thetaY(order);
                std::vector<double> thetaZ(order);
                
                // Calculate B-spline coefficients for each dimension
                std::vector<double> coefficients(order);
                
                // X dimension B-spline
                computeBSplineCoefficients(gridFractions[0], order, coefficients);
                
                // Only output B-spline coefficients in debug mode
                for (int i = 0; i < order; i++) {
                    thetaX[i] = coefficients[i];
                }
                
                if (platform::is_debug_mode()) {
                    double xMax = 0.0, xMin = 0.0, xSum = 0.0;
                    for (int i = 0; i < order; i++) {
                        if (i == 0) {
                            xMax = xMin = thetaX[i];
                        } else {
                            xMax = std::max(xMax, thetaX[i]);
                            xMin = std::min(xMin, thetaX[i]);
                        }
                        xSum += thetaX[i];
                    }
                
                    platform::log(LogLevel::DEBUG, "X axis B-spline coefficients (gridFraction=", gridFractions[0], "):");
                    for (int i = 0; i < order; i++) {
                        platform::log(LogLevel::DEBUG, "theta_x[", i, "] = ", thetaX[i]);
                    }
                    platform::log(LogLevel::DEBUG, "X weight range: [", xMin, ", ", xMax, "], sum: ", xSum);
                }
                
                // Y dimension B-spline
                computeBSplineCoefficients(gridFractions[1], order, coefficients);
                
                // Only output B-spline coefficients in debug mode
                for (int i = 0; i < order; i++) {
                    thetaY[i] = coefficients[i];
                }
                
                if (platform::is_debug_mode()) {
                    double yMax = 0.0, yMin = 0.0, ySum = 0.0;
                    for (int i = 0; i < order; i++) {
                        if (i == 0) {
                            yMax = yMin = thetaY[i];
                        } else {
                            yMax = std::max(yMax, thetaY[i]);
                            yMin = std::min(yMin, thetaY[i]);
                        }
                        ySum += thetaY[i];
                    }
                
                    platform::log(LogLevel::DEBUG, "Y axis B-spline coefficients (gridFraction=", gridFractions[1], "):");
                    for (int i = 0; i < order; i++) {
                        platform::log(LogLevel::DEBUG, "theta_y[", i, "] = ", thetaY[i]);
                    }
                    platform::log(LogLevel::DEBUG, "Y weight range: [", yMin, ", ", yMax, "], sum: ", ySum);
                }
                
                // Z dimension B-spline
                computeBSplineCoefficients(gridFractions[2], order, coefficients);
                
                // Only output B-spline coefficients in debug mode
                for (int i = 0; i < order; i++) {
                    thetaZ[i] = coefficients[i];
                }
                
                if (platform::is_debug_mode()) {
                    double zMax = 0.0, zMin = 0.0, zSum = 0.0;
                    for (int i = 0; i < order; i++) {
                        if (i == 0) {
                            zMax = zMin = thetaZ[i];
                        } else {
                            zMax = std::max(zMax, thetaZ[i]);
                            zMin = std::min(zMin, thetaZ[i]);
                        }
                        zSum += thetaZ[i];
                    }
                
                    platform::log(LogLevel::DEBUG, "Z axis B-spline coefficients (gridFraction=", gridFractions[2], "):");
                    for (int i = 0; i < order; i++) {
                        platform::log(LogLevel::DEBUG, "theta_z[", i, "] = ", thetaZ[i]);
                    }
                    platform::log(LogLevel::DEBUG, "Z weight range: [", zMin, ", ", zMax, "], sum: ", zSum);
                    
                    // xSum and ySum are out of scope, so use these variables
                    platform::log(LogLevel::DEBUG, "Three-dimensional weight product sum theoretical value: approximately 1.0");
                    
                    // Output nearest grid point and potential value (for testing)
                    platform::log(LogLevel::DEBUG, "Atom ", atom_index, " nearest grid point information:");
                    platform::log(LogLevel::DEBUG, "Grid index: (", gridIndices[0], ",", gridIndices[1], ",", gridIndices[2], ")");
                    
                    // Output potential value of the point and its surrounding grid points
                    platform::log(LogLevel::DEBUG, "Nearest grid point potential values:");
                    for (int dx = -1; dx <= 1; dx++) {
                        for (int dy = -1; dy <= 1; dy++) {
                            for (int dz = -1; dz <= 1; dz++) {
                                int xi = (gridIndices[0] + dx + nx) % nx;
                                int yi = (gridIndices[1] + dy + ny) % ny;
                                int zi = (gridIndices[2] + dz + nz) % nz;
                                int index = xi * ny * nz + yi * nz + zi;
                                
                                double pot_val = pgp_params.potentialGrid[index].real();
                                
                                if (dx == 0 && dy == 0 && dz == 0) {
                                    platform::log(LogLevel::DEBUG, "→ Center point (", xi, ",", yi, ",", zi, "): ", pot_val);
                                } else if (std::abs(pot_val) > 1e-6) {
                                    // Only output non-zero surrounding points
                                    platform::log(LogLevel::DEBUG, "Point (", xi, ",", yi, ",", zi, "): ", pot_val);
                                }
                            }
                        }
                    }
                }
                
                // Interpolate potential
                double potential = 0.0;
                
                // Loop through all B-spline support points
                for (int ix = 0; ix < order; ix++) {
                    int xindex = (gridIndices[0] + ix) % nx;
                    
                    for (int iy = 0; iy < order; iy++) {
                        int yindex = (gridIndices[1] + iy) % ny;
                        
                        for (int iz = 0; iz < order; iz++) {
                            int zindex = (gridIndices[2] + iz) % nz;
                            
                            // Calculate three-dimensional grid index
                            int index = xindex * ny * nz + yindex * nz + zindex;
                            
                            // Use B-spline weights to accumulate potential
                            double grid_value = pgp_params.potentialGrid[index].real();
                            double weight = thetaX[ix] * thetaY[iy] * thetaZ[iz];
                            potential += grid_value * weight;
                            
                            // Only output detailed information about important grid points in debug mode
                            if (platform::is_debug_mode() && std::abs(grid_value) > 1e-6 && ix < 2 && iy < 2 && iz < 2) {
                                platform::log(LogLevel::DEBUG, "Grid point (", xindex, ",", yindex, ",", zindex, ") potential=", grid_value, 
                                            ", weight=", weight, " (Each dimension weights: ", thetaX[ix], ",", thetaY[iy], ",", thetaZ[iz], ")");
                            }
                        }
                    }
                }
                
                // Only output weight sum in weight interpolation in debug mode
                if (platform::is_debug_mode()) {
                    double totalWeight = 0.0;
                    for (int ix = 0; ix < order; ix++) {
                        for (int iy = 0; iy < order; iy++) {
                            for (int iz = 0; iz < order; iz++) {
                                totalWeight += thetaX[ix] * thetaY[iy] * thetaZ[iz];
                            }
                        }
                    }
                    platform::log(LogLevel::DEBUG, "B-spline weight sum: ", totalWeight);
                }
                
                // Accumulate energy (potential * charge)
                double atom_energy = potential * atom.charge;
                raw_energy += atom_energy;
                
                if (platform::is_debug_mode()) {
                    platform::log(LogLevel::DEBUG, "Atom potential: ", potential, ", Atom energy contribution: ", atom_energy);
                }
            }
        }
    } else {
        // Normal processing for movement residues
        // ...
    }
    
    // Calculate energy summary
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Summing atomic contributions for total energy");
    }
    
    // Potential has been halved in precomputation, now multiply by 2x factor to calculate energy
    // According to pgp.md, energy should be 2 * Σ(q_i * φ(r_i))
    energy = 2.0 * raw_energy;
    
    // Only output energy calculation details in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Energy calculation details:");
        platform::log(LogLevel::DEBUG, "Raw energy (raw): ", raw_energy, " kJ/mol");
        platform::log(LogLevel::DEBUG, "Applying 2x factor: × 2.0");
        platform::log(LogLevel::DEBUG, "Final energy: ", energy, " kJ/mol");
        
        platform::log(LogLevel::DEBUG, "Raw PGP energy: ", raw_energy, " kJ/mol");
        platform::log(LogLevel::DEBUG, "Final PGP energy (×2): ", energy, " kJ/mol");
        
        // Output final energy value and debugging information
        platform::log(LogLevel::DEBUG, "Final calculated PGP energy: ", energy, " kJ/mol");
    }
}

/**
 * @brief Calculate moving molecule energy through interpolation and return calculation result
 * 
 * This is a wrapper function for interpolateMoleculeEnergy, directly returning the calculated energy value
 * Convenient for Python calls and testing.
 * 
 * @param state System state, containing information about moving molecules and precomputed potential grid
 * @return Calculated energy value
 */
double calculateMoleculeEnergy(model::MCState& state) {
    double energy = 0.0;
    interpolateMoleculeEnergy(state, energy);
    
    // Only output debugging output in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Final calculated molecule energy: ", energy);
    }
    
    // If no fixed residues, potential grid may need to be recomputed
    if (std::abs(energy) < 1e-10) {
        // Check if it's because no fixed residues are causing the calculation problem
        int fixed_count = 0;
        if (platform::is_debug_mode() || std::abs(energy) < 1e-10) {
            // Only calculate fixed residue count when checking reasons
            for (int i = 0; i < state.activeResidueCount; ++i) {
                if (state.residues[i].fixed && state.residues[i].active) {
                    fixed_count++;
                }
            }
        }
        
        if (fixed_count == 0) {
            platform::log(LogLevel::WARNING, "No fixed residues found, energy near zero!");
            
            if (platform::is_debug_mode()) {
                // Precompute potential grid for all residues
                platform::log(LogLevel::DEBUG, "Attempting to precompute potential grid for all residues...");
            }
            
            precomputeGridPotential(state, false);
            
            // Recalculate energy
            interpolateMoleculeEnergy(state, energy);
        }
    }
    
    // Ensure return correct sign and magnitude of energy value
    return energy;
}

double computeMoleculeEnergyGlobal(model::MCState& state, const std::vector<int>& movementResidues, const std::vector<int>& nearbyResidues, int threadIndex) {
    // If parameters are not initialized, return 0
    if (!pgp_params.initialized) {
        platform::log(LogLevel::WARNING, "PGP parameters not initialized, returning 0 energy");
        return 0.0;
    }
    
    // Initialize total energy to 0
    double totalEnergy = 0.0;
    
    // Only output logs in debug mode
    if (platform::is_debug_mode()) {
        // Record thread index used, can be used for log tracking in multi-threaded calculations
        platform::log(LogLevel::DEBUG, "Using thread index: ", threadIndex, " calculating PGP energy");
        
        // Record nearby residue count, can be used in some algorithm variants for short-range energy correction
        platform::log(LogLevel::DEBUG, "Considering nearby residue count: ", nearbyResidues.size());
    }
    
    // For direct space (short-range) energy correction, nearby residues can be considered
    double directSpaceCorrection = 0.0;
    if (!nearbyResidues.empty()) {
        // Only output logs in debug mode
        if (platform::is_debug_mode()) {
            platform::log(LogLevel::DEBUG, "Calculating direct space correction for nearby residues");
        }
        
        // Here direct space correction calculation can be implemented
        // But current PGP implementation mainly focuses on precomputed grid potential part
        // If full direct space correction is needed, it should be implemented separately
    }
    
    // Process different residue energy calculations
    if (movementResidues.empty()) {
        // Only output logs in debug mode
        if (platform::is_debug_mode()) {
            platform::log(LogLevel::DEBUG, "Calculating energy for all residues");
        }
        
        // Find all non-fixed active residues
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (state.residues[i].active && !state.residues[i].fixed) {
                // Directly call calculateMoleculeEnergy, it will call interpolateMoleculeEnergy internally
                totalEnergy = calculateMoleculeEnergy(state);
                break;
            }
        }
    } else {
        // If specified movement residues, we need to modify state's movementResidues
        // First backup original movementResidues
        auto originalMovementResidues = state.movementResidues;
        
        // Clear and set new movementResidues
        state.movementResidues.clear();
        model::MCMovementResidueInfo info;
        info.startIndex = movementResidues[0];
        info.activeCount = movementResidues.size();
        state.movementResidues.push_back(info);
        
        // Calculate energy
        totalEnergy = calculateMoleculeEnergy(state);
        
        // Restore original movementResidues
        state.movementResidues = originalMovementResidues;
    }
    
    // Add direct space correction (if any)
    totalEnergy += directSpaceCorrection;
    
    // Only output energy calculation result in debug mode
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "PGP energy calculation result: ", totalEnergy, 
                    " (Thread: ", threadIndex, ", Considering nearby residues: ", !nearbyResidues.empty(), ")");
    }
    
    return totalEnergy;
}

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
void computeRealSpacePGP(model::MCState& state, bool movement_only, bool store_in_residues) {
    const auto& box = state.info.box;
    auto& atoms = state.atoms;
    auto& residues = state.residues;
    const float cutoff2 = pgp_params.cutoff * pgp_params.cutoff;

    // Reset electrostatic energy
    for(auto& residue : residues) {
        if(residue.active) {
            residue.energy_elec = 0.0f;
        }
    }
    
    // Real-space total energy
    double real_space_total = 0.0;
    
    // Add debug information
    int debug_count = 0;
    const int max_debug_pairs = 5;

    // Loop over all residue pairs
    for(int r1 = 0; r1 < state.activeResidueCount; r1++) {
        if(!residues[r1].active) continue;
        if(movement_only) {
            bool in_movement = false;
            for(const auto& movementInfo : state.movementResidues) {
                if(r1 >= movementInfo.startIndex && 
                   r1 < movementInfo.startIndex + movementInfo.activeCount) {
                    in_movement = true;
                    break;
                }
            }
            if(!in_movement) continue;
        }
        
        for(int r2 = r1 + 1; r2 < state.activeResidueCount; r2++) {
            if(!residues[r2].active) continue;
            
            // 如果两个残基都是固定的，那么它们之间的相互作用已经包含在预计算的网格势能中
            if(residues[r1].fixed && residues[r2].fixed) continue;
            
            // Loop over atoms in each residue
            for(int i = residues[r1].atomStart; 
                i < residues[r1].atomStart + residues[r1].atomCount; i++) {
                // Ensure atom index is valid
                if(i >= state.activeAtomCount) continue;
                
                for(int j = residues[r2].atomStart; 
                    j < residues[r2].atomStart + residues[r2].atomCount; j++) {
                    // Ensure atom index is valid
                    if(j >= state.activeAtomCount) continue;
                    
                    float dx = atoms[i].x - atoms[j].x;
                    float dy = atoms[i].y - atoms[j].y;
                    float dz = atoms[i].z - atoms[j].z;
                    
                    // Apply periodic boundary conditions
                    dx -= box[0] * round(dx / box[0]);
                    dy -= box[1] * round(dy / box[1]);
                    dz -= box[2] * round(dz / box[2]);
                    
                    float r2 = dx*dx + dy*dy + dz*dz;
                    
                    // Skip pairs beyond cutoff
                    if(r2 > cutoff2) continue;
                    
                    // Compute energy
                    float r = sqrt(r2);
                    float qi = atoms[i].charge;
                    float qj = atoms[j].charge;
                    
                    // Skip neutral atoms
                    if(std::abs(qi) < 1e-6 || std::abs(qj) < 1e-6) continue;
                    
                    // Calculate real space contribution for PGP - only erfc part
                    double term = pgp_params.erfcApprox(r);
                    double pair_energy = qi * qj * term / r;
                    
                    // Print debug information
                    if (platform::is_debug_mode() && debug_count < max_debug_pairs) {
                        platform::log(LogLevel::DEBUG, 
                            "Debug energyPGP: Atom pair (", i, ",", j, "): ",
                            "r = ", r, " nm, ",
                            "q1*q2 = ", qi * qj, ", ",
                            "erfc term = ", term, ", ",
                            "energy = ", pair_energy,
                            ", with COULOMB = ", COULOMB * pair_energy, " kJ/mol");
                        debug_count++;
                    }

                    // Accumulate to total energy
                    real_space_total += pair_energy;
                    
                    // Decide how to store energy based on parameters
                    if (store_in_residues) {
                        // Each residue gets half of the pair interaction energy
                        residues[r1].energy_elec += pair_energy / 2.0f;
                        residues[r2].energy_elec += pair_energy / 2.0f;
                    }
                }
            }
        }
    }
    
    // Store total real-space energy (not yet multiplied by COULOMB)
    state.ewald_energy.real_space = real_space_total;
}

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
double computeSelfEnergyPGP(model::MCState& state, bool movement_only) {
    double self_energy = 0.0;
    double sum_q2 = 0.0;
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Computing self energy for PGP method");
        platform::log(LogLevel::DEBUG, "Movement only: ", movement_only);
    }
    
    // Sum up squares of charges
    if (movement_only) {
        // Only include moving residues
        for (const auto& movementInfo : state.movementResidues) {
            for (int i = movementInfo.startIndex; 
                 i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                if (!state.residues[i].active) continue;
                
                // Sum q^2 for all atoms in this movement residue
                for (int j = 0; j < state.residues[i].atomCount; j++) {
                    int atomIdx = state.residues[i].atomStart + j;
                    double q = state.atoms[atomIdx].charge;
                    sum_q2 += q * q;
                }
            }
        }
    } 
    else {
        // Sum q^2 for all atoms in the system
        for(int i = 0; i < state.activeAtomCount; i++) {
            double charge = state.atoms[i].charge;
            double q2 = charge * charge;
            sum_q2 += q2;
        }
    }
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Sum of q² = ", sum_q2);
    }
    
    // Self-energy formula: -ONE_4PI_EPS0 * alpha / sqrt(PI) * sum_q2
    double prefactor = -COULOMB * pgp_params.alpha / sqrt(M_PI);
    self_energy = prefactor * sum_q2;
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Self energy prefactor = ", prefactor, 
                     ", resulting self energy = ", self_energy);
    }
    
    return self_energy;
}

/**
 * @brief Use PGP method to calculate system energy
 * 
 * This function provides a complete energy calculation for the entire system using the
 * PGP-PME method, including grid potential interpolation, real-space electrostatics,
 * self energy correction, and Lennard-Jones interactions.
 * 
 * @param state MC state
 */
void computeSystemEnergyPGP(model::MCState& state) {
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized. Call setPGPParameters() first.");
    }
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Computing total system energy using PGP method");
    }
    
    // 1. 计算网格势能插值部分（通过预计算的势能场）
    double grid_energy = 0.0;
    interpolateMoleculeEnergy(state, grid_energy);
    
    // 2. 计算实空间部分（处理短程电荷相互作用）
    computeRealSpacePGP(state, false, true);
    
    // 3. 计算自能修正项
    state.ewald_energy.self = computeSelfEnergyPGP(state, false);
    
    // 4. 使用直接截断法计算LJ相互作用
    computeSystemVdwEnergyCutoff(state);
    
    // 将实空间能量乘以COULOMB常数，网格能量和自能已经包含该常数
    state.ewald_energy.real_space *= COULOMB;
    
    // 计算总的LJ能量
    double vdw_total = 0.0;
    for (const auto& residue : state.residues) {
        if (residue.active) {
            vdw_total += residue.energy_vdw;
        }
    }
    
    // 计算总能量：网格势能 + 实空间静电能 + 自能修正 + LJ能量
    state.ewald_energy.reciprocal = grid_energy;  // 将网格势能存储在reciprocal字段中
    state.ewald_energy.total = grid_energy + state.ewald_energy.real_space + 
                             state.ewald_energy.self + vdw_total;
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "PGP system energy components: ");
        platform::log(LogLevel::DEBUG, "  Grid energy = ", grid_energy);
        platform::log(LogLevel::DEBUG, "  Real space = ", state.ewald_energy.real_space);
        platform::log(LogLevel::DEBUG, "  Self energy = ", state.ewald_energy.self);
        platform::log(LogLevel::DEBUG, "  VDW energy = ", vdw_total);
        platform::log(LogLevel::DEBUG, "  Total energy = ", state.ewald_energy.total);
    }
}

/**
 * @brief Use PGP method to calculate energy of moving residues
 * 
 * This function calculates the energy for just the moving residues using the PGP-PME method,
 * which is useful for Monte Carlo move acceptance/rejection decisions.
 * 
 * @param state MC state
 */
void computeMovementEnergyPGP(model::MCState& state) {
    if (!pgp_params.initialized) {
        throw std::runtime_error("PGP parameters not initialized. Call setPGPParameters() first.");
    }
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "Computing movement residue energy using PGP method");
    }
    
    // 1. 计算网格势能插值部分（通过预计算的势能场）
    double grid_energy = 0.0;
    interpolateMoleculeEnergy(state, grid_energy);
    
    // 2. 计算实空间部分（处理短程电荷相互作用，仅计算移动残基）
    computeRealSpacePGP(state, true, true);
    
    // 3. 计算自能修正项（仅计算移动残基）
    state.ewald_energy.self = computeSelfEnergyPGP(state, true);
    
    // 4. 使用直接截断法计算LJ相互作用
    computeSystemVdwEnergyCutoff(state);
    
    // 将实空间能量乘以COULOMB常数
    state.ewald_energy.real_space *= COULOMB;
    
    // 只累计移动残基的LJ能量
    double vdw_total = 0.0;
    for (const auto& movementInfo : state.movementResidues) {
        for (int i = movementInfo.startIndex;
             i < movementInfo.startIndex + movementInfo.activeCount; i++) {
            if (state.residues[i].active) {
                vdw_total += state.residues[i].energy_vdw;
            }
        }
    }
    
    // 计算总能量
    state.ewald_energy.reciprocal = grid_energy;  // 将网格势能存储在reciprocal字段中
    state.ewald_energy.total = grid_energy + state.ewald_energy.real_space + 
                             state.ewald_energy.self + vdw_total;
    
    if (platform::is_debug_mode()) {
        platform::log(LogLevel::DEBUG, "PGP movement energy components: ");
        platform::log(LogLevel::DEBUG, "  Grid energy = ", grid_energy);
        platform::log(LogLevel::DEBUG, "  Real space = ", state.ewald_energy.real_space);
        platform::log(LogLevel::DEBUG, "  Self energy = ", state.ewald_energy.self);
        platform::log(LogLevel::DEBUG, "  VDW energy = ", vdw_total);
        platform::log(LogLevel::DEBUG, "  Total energy = ", state.ewald_energy.total);
    }
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 
