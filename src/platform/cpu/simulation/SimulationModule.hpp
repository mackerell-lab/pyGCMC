#pragma once

/**
 * @file SimulationModule.hpp
 * @brief Internal module interface for simulation components
 * 
 * This file provides internal interfaces and utilities for simulation modules,
 * similar to EnergyModule.hpp and MovementModule.hpp
 */

// Include all simulation sub-modules
#include "core/SimulationCore.hpp"
#include "core/SimulationRunner.hpp"
#include "setup/SystemInitializer.hpp"
#include "stats/StatisticsTracker.hpp"
#include "io/SimulationIO.hpp"
#include "gcmc/GCMCController.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

// Forward declarations for internal use
namespace core {
    class SimulationCore;
    class SimulationRunner;
}

namespace setup {
    class SystemInitializer;
}

namespace stats {
    class StatisticsTracker;
}

namespace io {
    class SimulationIO;
}

namespace gcmc {
    class GCMCController;
}

namespace impl {
    class GCMCSimulation;
    class GCMCSimulationModular;
}

// Internal utility functions
namespace internal {
    
    /**
     * @brief Create default simulation configuration
     */
    inline SimulationCore::Config createDefaultCoreConfig() {
        SimulationCore::Config config;
        config.numSteps = 100000;
        config.equilibrationSteps = 10000;
        config.collectStatistics = true;
        config.statisticsInterval = 100;
        config.energyMethod = EnergyMethod::DIRECT;
        config.cutoff = 12.0;
        config.randomSeed = -1;
        return config;
    }
    
    /**
     * @brief Create default I/O configuration
     */
    inline io::SimulationIO::Config createDefaultIOConfig() {
        io::SimulationIO::Config config;
        config.outputPrefix = "gcmc";
        config.trajectoryFormat = "pdb";
        config.compressTrajectory = false;
        config.writeEnergies = true;
        config.writeDensities = true;
        config.precision = 6;
        config.appendMode = false;
        return config;
    }
    
} // namespace internal

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc