#pragma once

/**
 * @file MovementModule.hpp
 * @brief Movement Module - Unified Entry Point for All GCMC Movement Functionality
 *
 * This is the ONLY header file you need to include to access all GCMC movement
 * capabilities in the simulation framework. The module provides comprehensive
 * GCMC functionality with all components implemented in C++ for maximum performance.
 *
 * **Module Organization:**
 *
 * **gcmc/ directory** - Complete GCMC Implementation
 * - GCMCModule.hpp: Main GCMC module with all functionality
 * - GCMCEngine.hpp: Core engine for move execution
 * - MoveSelector.hpp: Intelligent move selection
 * - BiasCalculator.hpp: Advanced biasing calculations
 * - GCMCStatistics.hpp: Comprehensive statistics
 * - AcceptanceCalculator.hpp: Acceptance criteria
 *
 * **reservoir/ directory** - Fragment Management
 * - FragmentReservoir.hpp: Fragment templates and instances
 *
 * **bias/ directory** - Advanced Biasing Techniques
 * - CavityBias.hpp: Cavity detection and biased insertion
 * - ConfigBias.hpp: Configurational bias (CBMC)
 *
 * **pool/ directory** - Memory Management
 * - ActivePool.hpp: Pre-allocated memory pool with ghost recycling
 *
 * **common/ directory** - Common Types and Utilities
 * - Vector3.hpp: 3D vector operations
 * - Quaternion.hpp: Rotation representations
 * - MovementUtils.hpp: Utility functions
 *
 * **Usage Example:**
 * @code
 * #include "platform/cpu/movement/MovementModule.hpp"
 * using namespace pygcmc::platform::cpu::movement;
 *
 * // Method 1: Use comprehensive GCMC module
 * gcmc::GCMCModule::Config config;
 * config.temperature = 300.0;
 * config.useCavityBias = true;
 * auto gcmcModule = gcmc::createGCMCModule(config);
 * gcmcModule->initialize(state);
 * gcmcModule->runProduction(100000);
 *
 * // Method 2: Use individual components
 * FragmentReservoir reservoir;
 * reservoir.addTemplate(waterTemplate);
 * GCMCEngine engine;
 * engine.initialize(&state, &reservoir);
 * auto result = engine.attemptInsertion(0);
 * @endcode
 */

// Core GCMC implementation - ALL IN C++
#include "gcmc/GCMCModule.hpp"
#include "gcmc/GCMCEngine.hpp"
#include "gcmc/GCMCMoveSelector.hpp"
#include "gcmc/GCMCBias.hpp"
#include "gcmc/GCMCStats.hpp"
#include "gcmc/GCMCAcceptance.hpp"

// Fragment reservoir system
#include "reservoir/FragmentReservoir.hpp"

// Advanced biasing techniques
#include "bias/CavityBias.hpp"
#include "bias/ConfigBias.hpp"

// Memory pool management
#include "pool/ActivePool.hpp"

// Common utilities - use model definitions instead
// Vector3 and Quaternion are defined in model/montecarlo/MCStructures.hpp

// Legacy movement components (for backward compatibility)
#include "core/MovementMain.hpp"
#include "common/MovementParams.hpp"
#include "common/MovementResult.hpp"
#include "common/MovementStatistics.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

/**
 * @brief Quick access to GCMC functionality
 *
 * These functions provide simple interfaces to run common GCMC simulations
 * without needing to configure all parameters manually.
 */
namespace GCMC {

    /**
     * @brief Run a standard GCMC simulation
     * @param state System state
     * @param fragments Fragment templates to insert/delete
     * @param config GCMC configuration
     * @param steps Number of steps to run
     * @return Final statistics
     */
    inline gcmc::GCMCStats runSimulation(
        model::montecarlo::MCState& state,
        const std::vector<FragmentTemplate>& fragments,
        const gcmc::GCMCModule::Config& config,
        int steps) {

        auto module = gcmc::createGCMCModule(config);
        module->initialize(state);

        for (const auto& fragment : fragments) {
            module->addFragmentType(fragment);
        }

        module->runEquilibration(config.equilibrationSteps);
        module->runProduction(steps);

        return module->getStatistics();
    }

    /**
     * @brief Create a water GCMC simulation
     */
    inline std::unique_ptr<gcmc::GCMCModule> createWaterGCMC(
        double temperature = 300.0,
        double chemicalPotential = -15.7) {

        gcmc::GCMCModule::Config config;
        config.temperature = temperature;
        config.useCavityBias = true;
        config.useConfigBias = false;

        auto module = gcmc::createGCMCModule(config);

        // Add water template
        FragmentTemplate water;
        water.name = "WAT";
        water.chemicalPotential = chemicalPotential;
        water.calculateActivity(temperature);
        // Note: Atoms should be added to water template

        return module;
    }

    /**
     * @brief Create a gas adsorption simulation
     */
    inline std::unique_ptr<gcmc::GCMCModule> createGasAdsorption(
        const std::string& gasType,
        double temperature = 298.0,
        double pressure = 1.0) {

        gcmc::GCMCModule::Config config;
        config.temperature = temperature;
        config.useCavityBias = true;

        auto module = gcmc::createGCMCModule(config);

        // Calculate chemical potential from pressure
        double mu = 8.314e-3 * temperature * std::log(pressure * 1e5);

        FragmentTemplate gas;
        gas.name = gasType;
        gas.chemicalPotential = mu;
        gas.calculateActivity(temperature);

        return module;
    }
}

/**
 * @brief Fragment template utilities
 */
namespace Templates {

    /**
     * @brief Create TIP3P water template
     */
    inline FragmentTemplate createTIP3P() {
        FragmentTemplate water;
        water.name = "WAT";

        // Add atoms
        model::montecarlo::MCAtom o;
        o.name = "O";
        o.type = 0;  // Type index for oxygen
        o.position = model::montecarlo::Vector3(0.0, 0.0, 0.0);
        o.charge = -0.834;
        o.mass = 15.999;
        water.atoms.push_back(o);

        model::montecarlo::MCAtom h1;
        h1.name = "H1";
        h1.type = 1;  // Type index for hydrogen
        h1.position = model::montecarlo::Vector3(0.0957, 0.0, 0.0);
        h1.charge = 0.417;
        h1.mass = 1.008;
        water.atoms.push_back(h1);

        model::montecarlo::MCAtom h2;
        h2.name = "H2";
        h2.type = 1;  // Type index for hydrogen
        h2.position = model::montecarlo::Vector3(-0.024, 0.0927, 0.0);
        h2.charge = 0.417;
        h2.mass = 1.008;
        water.atoms.push_back(h2);

        return water;
    }

    /**
     * @brief Create methane template
     */
    inline FragmentTemplate createMethane() {
        FragmentTemplate methane;
        methane.name = "CH4";

        model::montecarlo::MCAtom c;
        c.name = "C";
        c.type = 4;  // Type index for methane carbon
        c.position = model::montecarlo::Vector3(0.0, 0.0, 0.0);
        c.charge = -0.24;
        c.mass = 12.011;
        methane.atoms.push_back(c);

        // Add hydrogens in tetrahedral geometry
        // C-H bond length: 0.109 nm (not implemented yet)

        return methane;
    }

    /**
     * @brief Create CO2 template
     */
    inline FragmentTemplate createCO2() {
        FragmentTemplate co2;
        co2.name = "CO2";

        model::montecarlo::MCAtom c;
        c.name = "C";
        c.type = 2;  // Type index for CO2 carbon
        c.position = model::montecarlo::Vector3(0.0, 0.0, 0.0);
        c.charge = 0.70;
        c.mass = 12.011;
        co2.atoms.push_back(c);

        model::montecarlo::MCAtom o1;
        o1.name = "O1";
        o1.type = 3;  // Type index for CO2 oxygen
        o1.position = model::montecarlo::Vector3(-0.116, 0.0, 0.0);
        o1.charge = -0.35;
        o1.mass = 15.999;
        co2.atoms.push_back(o1);

        model::montecarlo::MCAtom o2;
        o2.name = "O2";
        o2.type = 3;  // Type index for CO2 oxygen
        o2.position = model::montecarlo::Vector3(0.116, 0.0, 0.0);
        o2.charge = -0.35;
        o2.mass = 15.999;
        co2.atoms.push_back(o2);

        return co2;
    }
}

// Export main classes for convenience
using GCMCModule = gcmc::GCMCModule;
using GCMCEngine = gcmc::GCMCEngine;
using GCMCMoveSelector = gcmc::GCMCMoveSelector;
using GCMCBias = gcmc::GCMCBias;
using GCMCStats = gcmc::GCMCStats;
using GCMCAcceptance = gcmc::GCMCAcceptance;

// Legacy compatibility - MovementModule is defined in core/MovementMain.hpp
// Remove self-referential alias

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc
