// MultiInsertionCBMC.hpp
// Multi-insertion CBMC implementation for parallel molecular insertions

#pragma once

#include "../../../../model/montecarlo/MCMain.hpp"
#include "../common/MovementParams.hpp"
#include <vector>
#include <array>
#include <random>
#include <memory>

// Forward declare CavityManager
namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
class CavityManager;
}}}}

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

// Configuration for multi-insertion
struct MultiInsertionConfig {
    int numTrialsPerRegion = 10;      // CBMC trials per region
    int maxParallelInsertions = 100;  // Max regions to attempt
    double minSeparation = 1.5;       // Min distance between regions (nm)
    bool useCavityBias = false;       // Enable cavity bias integration
    bool useGPUBatch = false;         // GPU acceleration (phase 3)
    double chemicalPotential = -15.7; // Chemical potential in kJ/mol
    
    // Parameters for fine control
    double displacementFraction = 0.5; // Fraction of region.radius used for sampling (0..1]
    bool useRegionVolume = true;       // Use region sampling volume for Veff; otherwise box volume
    
    // Independence control parameters (new)
    double cutoffNm = 1.2;            // Energy cutoff in nm
    double moleculeExtentNm = 0.15;   // Maximum molecular radius in nm (e.g., 0.15 for water)
    bool enforceIndependence = true;  // When true, enforce minimum separation for independence
    bool recomputeAfterAccept = false;// Fallback: recompute energies after each accept (expensive)
};

// Insertion region data
struct InsertionRegion {
    std::array<double, 3> center;     // Region center (nm)
    double radius;                     // Region radius (nm)
    std::array<int, 3> gridIndex;     // Grid position
    
    // Trial configurations
    std::vector<std::vector<model::montecarlo::MCAtom>> trialConfigs;
    std::vector<double> trialEnergies;
    
    // Selection results
    int selectedConfig = -1;
    double rosenbluthWeight = 0.0;
    bool accepted = false;
};

class MultiInsertionCBMC {
public:
    explicit MultiInsertionCBMC(const MultiInsertionConfig& config,
                                CavityManager* cavityManager = nullptr);
    ~MultiInsertionCBMC() = default;
    
    // Main interface
    std::pair<int, std::vector<InsertionRegion>> performMultiInsertion(
        model::montecarlo::MCState& state,
        int moleculeType,
        const MovementParams& params);
    
    // RNG control for reproducibility
    void setSeed(uint64_t seed);
    
    // Statistics
    double getAcceptanceRate() const;
    double getParallelEfficiency() const;
    void resetStatistics();
    
private:
    // Configuration
    MultiInsertionConfig config_;
    
    // Random number generation
    std::mt19937 rng_;
    std::uniform_real_distribution<> uniform_;
    
    // Cavity bias support
    CavityManager* cavityManager_;  // Non-owning pointer
    
    // Statistics tracking
    struct Statistics {
        int totalAttempts = 0;
        int totalAccepts = 0;
        std::vector<int> parallelAttempts;
        std::vector<int> parallelAccepts;
    };
    Statistics stats_;
    
    // Internal methods
    std::vector<InsertionRegion> divideBoxIntoRegions(const model::montecarlo::MCState& state);
    std::vector<InsertionRegion> selectNonAdjacentRegions(
        const std::vector<InsertionRegion>& allRegions,
        int numToSelect);
    
    void generateTrialConfigurations(
        std::vector<InsertionRegion>& regions,
        int moleculeType,
        const model::montecarlo::MCState& state);
    
    void batchCalculateEnergies(
        std::vector<InsertionRegion>& regions,
        model::montecarlo::MCState& state);
    
    void selectOptimalConfigurations(
        std::vector<InsertionRegion>& regions,
        double beta);
    
    std::vector<InsertionRegion> acceptInsertions(
        std::vector<InsertionRegion>& regions,
        const model::montecarlo::MCState& state,
        double beta,
        double chemicalPotential);
    
    // Utility functions
    std::vector<model::montecarlo::MCAtom> createMolecule(
        int moleculeType,
        const std::array<double, 3>& position);
    
    std::array<double, 4> generateRandomQuaternion();
    
    void rotateWithQuaternion(
        std::vector<model::montecarlo::MCAtom>& atoms,
        const std::array<double, 4>& quaternion);
    
    double calculateRosenbluthWeight(
        const std::vector<double>& energies,
        double beta);
    
    int selectByGumbelMax(const std::vector<double>& logWeights);
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc