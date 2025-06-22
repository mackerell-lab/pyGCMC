#pragma once
#include "model/ModelModule.hpp"

namespace pygcmc {
namespace system {
namespace montecarlo {

/**
 * @brief CHARMM-style switching function manager
 * 
 * Handles smooth energy cutoff using CHARMM's switching function
 * for proper energy conservation in MD/MC simulations.
 */
class MCSwitching {
public:
    MCSwitching() = default;
    ~MCSwitching() = default;

    // Disable copy operations
    MCSwitching(const MCSwitching&) = delete;
    MCSwitching& operator=(const MCSwitching&) = delete;

    // Enable move operations
    MCSwitching(MCSwitching&&) = default;
    MCSwitching& operator=(MCSwitching&&) = default;

    /**
     * @brief Set or disable CHARMM-style smooth switching function
     * 
     * @param state Monte Carlo state to modify
     * @param enable Whether to enable the switching function
     * @param r_on Inner cutoff radius (nm), distance at which attenuation begins
     * @param r_off Outer cutoff radius (nm), distance at which energy drops to zero
     * 
     * @throw std::runtime_error If parameters are invalid (r_on >= r_off or negative values)
     */
    void setSwitchingFunction(model::MCState& state, bool enable, float r_on = 1.0f, float r_off = 1.2f);

    /**
     * @brief Calculate the switching function value at the given distance
     * 
     * Uses CHARMM's switching function:
     * S(r) = [(r_off^2 - r^2)^2 * (r_off^2 + 2r^2 - 3r_on^2)] / (r_off^2 - r_on^2)^3
     * 
     * @param state Monte Carlo state containing switching parameters
     * @param r Distance (nm) at which to calculate the switching function value
     * @return Switching function value, ranging from [0,1]
     */
    float calculateSwitchingFunction(const model::MCState& state, float r) const;

    /**
     * @brief Get whether the switching function is currently enabled
     * 
     * @param state Monte Carlo state to check
     * @return Whether the switching function is enabled
     */
    bool isUsingSwitchingFunction(const model::MCState& state) const;

    /**
     * @brief Get inner cutoff radius
     * 
     * @param state Monte Carlo state to check
     * @return Inner cutoff radius (nm)
     */
    float getSwitchingROn(const model::MCState& state) const;

    /**
     * @brief Get outer cutoff radius
     * 
     * @param state Monte Carlo state to check
     * @return Outer cutoff radius (nm)
     */
    float getSwitchingROff(const model::MCState& state) const;

    /**
     * @brief Apply current switching function settings to external state object
     * 
     * @param sourceState Source state with switching parameters
     * @param targetState Target state to apply settings to
     */
    void applySwitchingToState(const model::MCState& sourceState, model::MCState& targetState) const;
};

} // namespace montecarlo
} // namespace system
} // namespace pygcmc 