// src/simulation/simulation.hpp

#pragma once

#include "../platform/platform.hpp"
#include "../model/montecarlo.hpp"
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
     * 设置PGP-PME (Precomputed Grid-Potential Particle Mesh Ewald)算法参数
     * 
     * 配置用于加速蒙特卡洛模拟中电荷相互作用计算的PGP-PME算法参数。
     * 该方法优化了传统PME方法，特别适用于MC模拟中移动部分与固定部分间的相互作用计算。
     * 
     * @param alpha Ewald分离参数
     * @param meshSize PME网格大小
     * @param pair_cutoff 配对相互作用截断距离
     * @param pairGridSize 预计算电势网格大小
     * @param splineOrder B样条插值阶数
     * @param tolerance 误差容限
     */
    void setPGPParameters(double alpha, const std::array<int, 3>& meshSize, double pair_cutoff,
                          const std::array<int, 3>& pairGridSize, int splineOrder, double tolerance);

private:
    std::unique_ptr<platform::IPlatform> platform_;
};

} // namespace simulation
} // namespace pygcmc