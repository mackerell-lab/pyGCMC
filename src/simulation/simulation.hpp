// src/simulation/simulation.hpp

#pragma once

#include "../platform/platform.hpp"
#include "../model/montecarlo.hpp"
#include <memory>
#include <sstream>
#include <iostream>

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
static bool verbose_ = false;
static LogLevel log_level_ = LogLevel::INFO;

// Logging functions
inline void set_verbose(bool verbose) { verbose_ = verbose; }
inline void set_log_level(LogLevel level) { log_level_ = level; }

template<typename... Args>
inline void log(LogLevel level, Args... args) {
    if (!verbose_ || level < log_level_) return;
    
    std::stringstream ss;
    (ss << ... << args);
    
    switch (level) {
        case LogLevel::DEBUG:
            std::cout << "[SIMULATION DEBUG] ";
            break;
        case LogLevel::INFO:
            std::cout << "[SIMULATION INFO] ";
            break;
        case LogLevel::WARNING:
            std::cout << "[SIMULATION WARNING] ";
            break;
        case LogLevel::ERROR:
            std::cout << "[SIMULATION ERROR] ";
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

    // Finalize the simulation and (optionally) download final data
    void finalize() {
        platform_->finalize();
    }

private:
    std::unique_ptr<platform::IPlatform> platform_;
};

} // namespace simulation
} // namespace pygcmc