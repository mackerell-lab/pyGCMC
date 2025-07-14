#pragma once

#include "model/ModelModule.hpp"
#include "platform/platform.hpp"
#include "../common/EnergyConstants.hpp"
#include "../common/EnergyUtils.hpp"
#include <array>
#include <vector>
#include <complex>
#include <mutex>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Independent PGP Parameters - No inheritance from PME
 * 
 * This structure contains all parameters required by PGP algorithm
 * without inheriting from PME, preventing shared state issues.
 */
struct PGPParamsIndependent {
    // Basic parameters (previously inherited from PME)
    double alpha{1.0};                 // Ewald separation parameter
    int meshSize[3]{64,64,64};         // Mesh size for FFT
    double tolerance{1e-5f};           // Error tolerance
    bool initialized{false};           // Initialization flag
    double cutoff{0.0};               // Real space cutoff
    double epsilon_r{1.0};            // Relative dielectric constant
    double box[3]{1.0, 1.0, 1.0};     // Box dimensions
    
    // Independent lookup tables (deep copied, not shared)
    std::vector<double> erfcTable;         // Own erfc lookup table
    std::vector<double> ewaldScaleTable;   // Own Ewald scaling table
    double ewaldDX;                        // Table step size
    double ewaldDXInv;                     // Inverse step size
    double erfcDXInv;                      // Inverse erfc step size
    
    // B-spline parameters (independent)
    int splineOrder{4};                    // B-spline order
    std::vector<double> bsplineModuli[3];  // Own B-spline moduli
    
    // Independent grid structures (no sharing with PME!)
    std::vector<std::complex<double>> pmeGrid;    // Own FFT grid
    std::vector<double> pmeCharge;                // Own charge grid
    
    // PGP-specific parameters
    double potential_cutoff;                       // Potential cutoff
    int potential_grid_size[3];                    // Potential grid dimensions
    double grid_spacing;                           // Grid spacing
    std::vector<std::complex<double>> potentialGrid;  // Precomputed potential
    
    // Debug flag
    bool debug_mode = true;
    
    /**
     * Initialize all internal data structures independently
     */
    void initializeIndependent() {
        // Initialize own lookup tables
        initializeOwnTables();
        
        // Initialize own B-splines
        initializeOwnBsplines();
        
        // Initialize own grids
        initializeOwnGrids();
        
        // Initialize potential grid
        initializePotentialGrid();
        
        initialized = true;
    }
    
private:
    void initializeOwnTables();
    void initializeOwnBsplines();
    void initializeOwnGrids();
    void initializePotentialGrid();
};

// New global PGP parameters - completely independent
extern PGPParamsIndependent pgp_params_independent;

// Independent initialization functions
void setPGPParametersIndependent(double alpha, const int meshSize[3], 
                                double potential_cutoff, 
                                const int potentialGridSize[3], 
                                int splineOrder, double tolerance);

void initializePGPParametersIndependent(double cutoff, const double box[3], 
                                       double alpha, const int meshSize[3], 
                                       double potentialCutoff,
                                       const int potentialGridSize[3],
                                       int splineOrder, double tolerance);

// Functions that use independent parameters
void computeSystemEnergyPGPIndependent(model::MCState& state);
void computeMovementEnergyPGPIndependent(model::MCState& state);

} // namespace cpu
} // namespace platform
} // namespace pygcmc