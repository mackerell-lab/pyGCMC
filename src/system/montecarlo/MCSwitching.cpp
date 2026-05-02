#include "MCSwitching.hpp"
#include <stdexcept>

namespace pygcmc {
namespace system {
namespace montecarlo {

void MCSwitching::setSwitchingFunction(model::MCState& state, bool enable, float r_on, float r_off) {
    // Parameter validation
    if (r_on >= r_off) {
        throw std::runtime_error("Invalid switching function parameters: r_on must be less than r_off");
    }
    if (r_on <= 0.0f || r_off <= 0.0f) {
        throw std::runtime_error("Invalid switching function parameters: radii must be positive");
    }

    // Set parameters in MCState
    state.info.use_switching = enable;
    state.info.r_on = r_on;
    state.info.r_off = r_off;
}

float MCSwitching::calculateSwitchingFunction(const model::MCState& state, float r) const {
    if (!state.info.use_switching || r <= state.info.r_on) {
        return 1.0f;  // No switching function applied when r <= r_on
    }
    if (r >= state.info.r_off) {
        return 0.0f;  // Energy is zero when r >= r_off
    }

    // Calculate CHARMM-style switching function
    // S(r) = [(r_off^2 - r^2)^2 * (r_off^2 + 2r^2 - 3r_on^2)] / (r_off^2 - r_on^2)^3
    float r2 = r * r;
    float ron2 = state.info.r_on * state.info.r_on;
    float roff2 = state.info.r_off * state.info.r_off;

    float numerator = (roff2 - r2) * (roff2 - r2) * (roff2 + 2.0f*r2 - 3.0f*ron2);
    float denominator = (roff2 - ron2) * (roff2 - ron2) * (roff2 - ron2);

    return numerator / denominator;
}

bool MCSwitching::isUsingSwitchingFunction(const model::MCState& state) const {
    return state.info.use_switching;
}

float MCSwitching::getSwitchingROn(const model::MCState& state) const {
    return state.info.r_on;
}

float MCSwitching::getSwitchingROff(const model::MCState& state) const {
    return state.info.r_off;
}

void MCSwitching::applySwitchingToState(const model::MCState& sourceState, model::MCState& targetState) const {
    targetState.info.use_switching = sourceState.info.use_switching;
    targetState.info.r_on = sourceState.info.r_on;
    targetState.info.r_off = sourceState.info.r_off;
}

} // namespace montecarlo
} // namespace system
} // namespace pygcmc
