#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_PARAMS_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_PARAMS_HPP

#include <cmath>
#include <cstdint>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

/**
 * Parameters for GCMC movement operations
 */
struct MovementParams {
    // Basic thermodynamic parameters
    double temperature = 298.15;           // Temperature in Kelvin
    double beta = 1.0 / (8.314e-3 * 298.15);  // 1/kT in mol/kJ
    double chemicalPotential = -15.7;      // Chemical potential in kJ/mol
    
    // Cavity bias parameters
    bool useCavityBias = true;             // Enable cavity bias for insertion
    double cavityGridSpacing = 0.2;        // Grid spacing in nm (was 2.0 Å)
    double probeRadius = 0.14;             // Probe radius for cavity detection in nm (was 1.4 Å)
    int cavityUpdateFrequency = 100;       // Update cavity grid every N moves
    
    // Configurational bias parameters
    bool useConfigBias = true;             // Enable configurational bias
    bool useConfigBiasForInsertion = false;// Enable CBMC for insertion (two-step method)
    int numConfigTrials = 10;              // Number of trial configurations (reduced for speed)
    bool includeTranslationInConfig = true;// Also vary position in config bias
    double configTranslationRange = 0.05;  // Translation range for config bias in nm (was 0.5 Å)
    
    // Movement limits
    double maxTranslation = 0.1;           // Maximum translation distance in nm (was 1.0 Å)
    double maxRotation = M_PI;             // Maximum rotation angle in radians
    
    // Numerical stability parameters
    bool useLogSpace = true;               // Use log-space calculations for stability
    double logSpaceMin = 1e-10;           // Minimum probability value
    double logSpaceMax = 1e10;            // Maximum probability value
    
    // Two-stage strategy parameters
    bool useTwoStageStrategy = false;      // Enable two-stage insertion strategy
    double pgpThreshold = -10.0;           // PGP energy threshold in kJ/mol
    int maxCandidatesStage2 = 16;          // Max candidates for stage 2 evaluation
    
    // Active pool parameters
    int maxAtoms = 100000;                 // Maximum atoms in active pool
    int maxResidues = 30000;               // Maximum residues in active pool
    double fragmentationThreshold = 0.3;   // Trigger compaction when fragmentation > threshold
    
    // System parameters
    double volumeNm3 = 0.0;                // System volume in nm^3 (computed from box)
    double idealGasConcentration = 0.0;    // Ideal gas concentration (computed)
    
    // Random number generator seed
    uint64_t seed = 0;                     // RNG seed (0 = use time-based seed)
    
    // Movement probabilities (should sum to 1.0)
    double insertionProbability = 0.25;
    double deletionProbability = 0.25;
    double translationProbability = 0.25;
    double rotationProbability = 0.25;
    
    // Multi-insertion CBMC parameters
    bool useMultiInsertionCBMC = false;    // Enable multi-insertion mode
    int maxParallelInsertions = 0;         // 0 disables, >0 enables
    double minRegionSeparationNm = 1.5;    // Region spacing (nm), must be >= cutoff
    
    // Update derived parameters
    void updateDerivedParameters() {
        beta = 1.0 / (8.314e-3 * temperature);  // Update beta when temperature changes
    }
    
    // Constructor with default values
    MovementParams() {
        updateDerivedParameters();
    }
    
    // Constructor with temperature
    explicit MovementParams(double temp) : temperature(temp) {
        updateDerivedParameters();
    }
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_PARAMS_HPP