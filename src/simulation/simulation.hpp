// src/simulation/simulation.hpp

#pragma once

#include "../platform/platform.hpp"
#include "../model/ModelModule.hpp"
#include <memory>
#include <sstream>
#include <iostream>
#include <vector>

namespace pygcmc {
namespace simulation {

// Log level enum
enum class LogLevel {
    DEBUG,
    INFO,
    WARNING,
    ERROR
};

// Static logging control
static bool verbose_ = false;  // Default to false for production use
static LogLevel log_level_ = LogLevel::WARNING;  // Default to WARNING level

// Logging functions
inline void set_verbose(bool verbose) { verbose_ = verbose; }
inline void set_log_level(LogLevel level) { log_level_ = level; }
inline void set_debug_mode(bool debug_mode) { platform::set_debug_mode(debug_mode); }

// Helper function to check if debug output is enabled
inline bool is_debug_enabled() { 
    return verbose_ && log_level_ <= LogLevel::DEBUG; 
}

template<typename... Args>
inline void log(LogLevel level, Args... args) {
    if (!verbose_ || level < log_level_) return;
    
    std::stringstream ss;
    (ss << ... << args);
    
    switch (level) {
        case LogLevel::DEBUG:
            std::cout << "[DEBUG] ";
            break;
        case LogLevel::INFO:
            std::cout << "[INFO] ";
            break;
        case LogLevel::WARNING:
            std::cout << "[WARN] ";
            break;
        case LogLevel::ERROR:
            std::cout << "[ERROR] ";
            break;
    }
    std::cout << ss.str() << std::endl;
}

class Simulation {
public:
    // Construct the simulation with a platform (e.g. CPU or CUDA)
    explicit Simulation(std::unique_ptr<platform::IPlatform> platform)
      : platform_(std::move(platform)) {}

    // Upload the initial state to the platform
    void initialize(const model::MCState& initialState) {
        platform_->initialize(initialState);
    }

    // Run the simulation for a number of steps by attempting moves
    void run(int steps) {
        for (int i = 0; i < steps; ++i) {
            // For simplicity, using "translate" as the move type
            platform_->attemptMove("translate");
        }
    }

    // Compute and return the total energy from the platform
    float computeEnergy() {
        return platform_->computeTotalEnergy();
    }

    // Compute naive nonbonded energy between active movement residues and all other active residues
    // Note: This function modifies the energy_vdw and energy_elec parameters of residues in the state
    static void computeNaiveNonbondedEnergy(model::MCState& state);

    /**
     * @brief Calculate nonbonded energies for movement residues only
     */
    static void computeMovementEnergy(model::MCState& state);

    /**
     * @brief Calculate nonbonded energies for movement residues only with distance cutoff
     */
    static void computeMovementEnergyCutoff(model::MCState& state);

    /**
     * @brief Calculate nonbonded energies for the full system
     */
    static void computeSystemEnergy(model::MCState& state);

    /**
     * @brief Calculate nonbonded energies for the full system with distance cutoff
     */
    static void computeSystemEnergyCutoff(model::MCState& state);

    /**
     * @brief Calculate nonbonded energies for the full system with periodic boundary conditions
     * 
     * This function calculates nonbonded interactions (VDW and electrostatic)
     * between all active residues without distance cutoff,
     * applying periodic boundary conditions using the minimum image convention.
     * 
     * @param state System state containing residues and force field parameters
     * @throws std::runtime_error if box dimensions are invalid for PBC calculation
     */
    static void computeSystemEnergyPBC(model::MCState& state);

    /**
     * @brief Calculate nonbonded energies for the full system with periodic boundary conditions and cutoff
     * 
     * This function calculates nonbonded interactions (VDW and electrostatic)
     * between all active residues within the specified cutoff distance,
     * applying periodic boundary conditions using the minimum image convention.
     * 
     * @param state System state containing residues and force field parameters
     * @throws std::runtime_error if box dimensions are invalid for PBC calculation
     */
    static void computeSystemEnergyPBCCutoff(model::MCState& state);

    static void setEnergyDebugOutput(bool enable);

    // Finalize the simulation and (optionally) download final data
    void finalize() {
        platform_->finalize();
    }

    // Ewald method interfaces
    static void setEwaldParameters(float alpha, const int kmax[3], float tolerance = 1e-5f);
    static void computeSystemEnergyEwald(model::MCState& state);
    static void computeMovementEnergyEwald(model::MCState& state);
    
    // PME method interfaces
    static void setPMEParameters(float alpha, const int meshSize[3], int splineOrder = 4, float tolerance = 1e-5f);
    static void computeSystemEnergyPME(model::MCState& state);
    static void computeMovementEnergyPME(model::MCState& state);
    
    // PGP method interfaces
    static void setPGPParameters(float alpha, const int meshSize[3], float pair_cutoff, 
                                  const int pairGridSize[3], int splineOrder = 4, float tolerance = 1e-5f);

    // Energy calculation methods exposed to Python
    static void computeSystemVdwEnergyCutoff(model::MCState& state);
    
    // CHARMM switching function related methods have been moved to MonteCarloSystem class
    // Commented out here to prevent users from using deprecated interfaces
    
    // Ewald methods
    static void initializeEwaldParameters(float cutoff, const float box[3], 
                                          float alpha = 0.0f, float tolerance = 1e-5f);
                                          
    // PME methods
    static void initializePMEParameters(float cutoff, const float box[3], 
                                       float alpha = 0.0f, const int* meshSize = nullptr,
                                       int splineOrder = 4, float tolerance = 1e-5f);
    
    // PGP methods
    static void initializePGPParameters(float cutoff, float pair_cutoff, const float box[3], 
                                         float alpha = 0.0f, const int* meshSize = nullptr,
                                         const int* pairGridSize = nullptr,
                                         int splineOrder = 4, float tolerance = 1e-5f);

    /**
     * Set PGP-PME (Precomputed Grid-Potential Particle Mesh Ewald) algorithm parameters
     * 
     * Configure parameters for PGP-PME algorithm which accelerates electrostatic interaction calculations in Monte Carlo simulations.
     * This method optimizes traditional PME method, specially designed for interaction calculations between moving and fixed parts in MC simulations.
     * 
     * @param alpha Ewald separation parameter
     * @param meshSize PME grid size
     * @param pair_cutoff Pair interaction cutoff distance
     * @param pairGridSize Precomputed potential grid size
     * @param splineOrder B-spline interpolation order
     * @param tolerance Error tolerance
     */
    void setPGPParameters(double alpha, const std::array<int, 3>& meshSize, double pair_cutoff,
                          const std::array<int, 3>& pairGridSize, int splineOrder, double tolerance);

    // PGP method - add the following functions
    
    /**
     * Reset PGP global state to fix memory corruption issues
     * 
     * This function clears all global PGP state and should be called
     * between test runs or when reinitializing PGP parameters
     */
    static void resetPGPState();
    
    /**
     * Precompute grid potential for fixed parts of the system
     * 
     * Precompute potential grid for PGP-PME algorithm, which is a key step for accelerating MC simulations
     * 
     * @param state System state
     * @param fixed_only Whether to process only fixed parts
     */
    static void precomputeGridPotential(model::MCState& state, bool fixed_only);
    
    /**
     * Calculate moving molecule energy through interpolation (output parameter version)
     * 
     * Interpolate energy of moving molecule from precomputed potential grid
     * 
     * @param state System state
     * @param energy Output parameter that stores the calculated energy
     */
    static void interpolateMoleculeEnergy(model::MCState& state, double& energy);
    
    /**
     * Calculate moving molecule energy through interpolation (return value version)
     * 
     * Interpolate energy of moving molecule from precomputed potential grid and return the result
     * 
     * @param state System state
     * @return The calculated energy value
     */
    static double calculateMoleculeEnergy(model::MCState& state);

    /**
     * @brief Use PGP method to calculate system energy
     * 
     * This function calculates energy for the entire system using the PGP-PME method.
     * 
     * @param state MC state
     */
    static void computeSystemEnergyPGP(model::MCState& state);

    /**
     * @brief Use PGP method to calculate energy of moving residues
     * 
     * This function calculates energy for just the moving residues using the PGP-PME method.
     * 
     * @param state MC state
     */
    static void computeMovementEnergyPGP(model::MCState& state);

    /**
     * @brief Use fixed PGP method to calculate system energy
     * 
     * This function uses the corrected PGP implementation that properly uses pgp_params
     * for erfc calculations and includes intra-residue interactions.
     * 
     * @param state MC state
     */
    static void computeSystemEnergyPGPFixed(model::MCState& state);

    /**
     * @brief Use fixed PGP method to calculate energy of moving residues
     * 
     * This function uses the corrected PGP implementation for moving residues.
     * 
     * @param state MC state
     */
    static void computeMovementEnergyPGPFixed(model::MCState& state);

    /**
     * @brief Clear PME engine state
     * 
     * This function clears all global PME state including parameters, grids,
     * FFT weights, and other cached data to prevent cross-test contamination.
     */
    static void clearPMEEngine();
    
    /**
     * @brief Compute complete system energy using PGP with all interactions
     * 
     * This version includes intramolecular LJ interactions similar to PME Complete,
     * providing a complete energy calculation for PGP method.
     * 
     * @param state MC state
     */
    static void computeSystemEnergyPGPComplete(model::MCState& state);
    
    /**
     * @brief Compute movement energy using PGP Complete
     * 
     * Calculates energy for movement residues including all LJ interactions.
     * 
     * @param state MC state
     */
    static void computeMovementEnergyPGPComplete(model::MCState& state);
    
    /**
     * @brief Corrected implementation of movement energy using PGP Complete
     * 
     * Properly calculates energy for movement residues:
     * - Grid interpolation for movement atoms
     * - Real space interactions with ALL atoms  
     * - LJ interactions with ALL atoms
     * - Self energy of movement atoms only
     * 
     * @param state MC state
     */
    static void computeMovementEnergyPGPCompleteCorrect(model::MCState& state);

    /**
     * @brief Compute complete system energy using PME with all interactions
     * 
     * This version includes intramolecular LJ interactions that are normally
     * excluded in the Fixed methods, providing energy values that match
     * reference implementations like OpenMM.
     * 
     * @param state MC state containing system information
     */
    static void computeSystemEnergyPMEComplete(model::MCState& state);

    /**
     * @brief Compute complete system energy using cutoff with all interactions
     * 
     * This version includes intramolecular LJ interactions for comparison
     * with PME complete method and reference implementations.
     * 
     * @param state MC state containing system information
     */
    static void computeSystemEnergyCutoffComplete(model::MCState& state);

    /**
     * @brief Get the total electrostatic and van der Waals energy components
     * 
     * This function sums the energy components from all active residues and returns
     * them as a pair. Note that the energies are already divided by 2 to account
     * for double counting in pairwise calculations.
     * 
     * @param state MC state containing system information
     * @return std::pair<double, double> A pair of (electrostatic_energy, vdw_energy)
     */
    static std::pair<double, double> getTotalEnergyComponents(const model::MCState& state);

private:
    std::unique_ptr<platform::IPlatform> platform_;
};

/**
 * @brief Simple GCMC interface for high-level usage
 * 
 * This provides a minimal API for running GCMC simulations without
 * exposing the complexity of the underlying implementation.
 */
class GCMCSimulation {
public:
    // Simple configuration
    struct Config {
        double temperature;      // K
        int equilibrationSteps;
        int productionSteps;
        double chemicalPotential; // kJ/mol (for single component)
        bool useCavityBias;
        bool verbose;
        
        // Constructor with default values
        Config()
            : temperature(300.0),
              equilibrationSteps(10000),
              productionSteps(100000),
              chemicalPotential(-15.7),
              useCavityBias(true),
              verbose(false) {}
    };
    
    // Constructor
    explicit GCMCSimulation(const Config& config = Config());
    ~GCMCSimulation();
    
    // Initialize with state
    void initialize(model::MCState& state);
    
    // Add water molecules (convenience method)
    void addWater();
    
    // Run simulation
    void run();
    void runSteps(int nSteps);
    
    // Get results
    struct Results {
        double averageMolecules;
        double averageEnergy;
        double acceptanceRate;
    };
    Results getResults() const;
    
private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

/**
 * @brief Quick GCMC functions for common use cases
 */
namespace GCMC {
    // Run water GCMC with default parameters
    void runWaterSimulation(
        model::MCState& state,
        double temperature = 300.0,
        double chemicalPotential = -15.7,
        int steps = 100000
    );
    
    // Calculate chemical potential from pressure (ideal gas)
    double pressureToChemicalPotential(
        double pressure,      // bar
        double temperature,   // K
        double molecularMass  // g/mol
    );
}

} // namespace simulation
} // namespace pygcmc