// MultiInsertionCBMC.hpp
// Multi-insertion CBMC implementation for parallel molecular insertions

#pragma once

#include "../../../../model/montecarlo/MCMain.hpp"
#include "../common/MovementParams.hpp"
#include <vector>
#include <array>
#include <random>
#include <memory>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

// Configuration for multi-insertion
struct MultiInsertionConfig {
    int numTrialsPerRegion = 10;      // CBMC trials per region
    int maxParallelInsertions = 100;  // Max regions to attempt
    double minSeparation = 15.0;      // Min distance between regions (Angstrom)
    bool useCavityBias = false;       // Enable cavity bias (phase 2)
    bool useGPUBatch = false;         // GPU acceleration (phase 3)
    double chemicalPotential = -15.7; // Chemical potential in kJ/mol
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
    explicit MultiInsertionCBMC(const MultiInsertionConfig& config);
    ~MultiInsertionCBMC() = default;
    
    // Main interface
    std::pair<int, std::vector<InsertionRegion>> performMultiInsertion(
        model::montecarlo::MCState& state,
        int moleculeType,
        const MovementParams& params);
    
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