#pragma once

/**
 * @file SimulationCore.hpp
 * @brief Core simulation engine for GCMC
 */

#include "../../../../model/montecarlo/MCMain.hpp"
#include "../../movement/gcmc/GCMCEngine.hpp"
#include "../../movement/MovementModule.hpp"
#include "../../energy/EnergyModule.hpp"
#include <memory>
#include <functional>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {
namespace core {

/**
 * @brief Core simulation engine that manages the main simulation loop
 */
class SimulationCore {
public:
    /**
     * @brief Callback types for simulation events
     */
    using StepCallback = std::function<void(int step, const model::montecarlo::MCState& state)>;
    using MoveCallback = std::function<void(const platform::cpu::movement::MovementResult& result)>;

    /**
     * @brief Configuration for the simulation core
     */
    struct Config {
        int numSteps = 100000;
        int equilibrationSteps = 10000;
        bool collectStatistics = true;
        int statisticsInterval = 100;
        EnergyMethod energyMethod = EnergyMethod::DIRECT;
        double cutoff = 12.0;
        int randomSeed = -1;

        Config() = default;
    };

    // Constructor and destructor
    explicit SimulationCore(const Config& config);
    ~SimulationCore();

    // Initialize with state and engine
    void initialize(model::montecarlo::MCState* state,
                   movement::gcmc::GCMCEngine* engine);

    // Run simulation
    bool run(int nSteps = 0);

    // Step-by-step execution
    bool performStep();
    bool performMove();

    // Control
    void stop();
    bool isRunning() const;
    void reset();

    // Callbacks
    void setStepCallback(StepCallback callback);
    void setMoveCallback(MoveCallback callback);

    // Statistics
    int getCurrentStep() const;
    int getTotalMoves() const;
    int getAcceptedMoves() const;
    double getAcceptanceRate() const;

    // Energy computation
    void computeSystemEnergy();
    void computeMovementEnergy();
    std::pair<double, double> getEnergyComponents() const;

    // Configuration
    const Config& getConfig() const { return config_; }
    void updateConfig(const Config& config);

private:
    Config config_;
    model::montecarlo::MCState* state_;
    movement::gcmc::GCMCEngine* engine_;

    bool initialized_;
    bool running_;
    bool stopRequested_;

    int currentStep_;
    int totalMoves_;
    int acceptedMoves_;

    StepCallback stepCallback_;
    MoveCallback moveCallback_;

    // Private methods
    void checkInitialized() const;
    void updateStatistics(const platform::cpu::movement::MovementResult& result);
    void applyEnergyMethod();
};

} // namespace core
} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc
