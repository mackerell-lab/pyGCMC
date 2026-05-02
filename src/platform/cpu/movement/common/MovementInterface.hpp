#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_INTERFACE_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_INTERFACE_HPP

#include <vector>
#include "MovementParams.hpp"
#include "MovementResult.hpp"

// Forward declaration
namespace pygcmc {
namespace model {
namespace montecarlo {
    class MCState;
}
}
}

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using model::montecarlo::MCState;

/**
 * Abstract interface for all movement operations
 */
class MovementInterface {
public:
    virtual ~MovementInterface() = default;

    // Core movement operations
    virtual MovementResult attemptInsertion(MCState& state, const MovementParams& params) = 0;
    virtual MovementResult attemptDeletion(MCState& state, const MovementParams& params) = 0;
    virtual MovementResult attemptTranslation(MCState& state, const MovementParams& params) = 0;
    virtual MovementResult attemptRotation(MCState& state, const MovementParams& params) = 0;

    // Optional: batch operations for efficiency
    virtual std::vector<MovementResult> attemptBatchInsertions(
        MCState& state,
        const MovementParams& params,
        int numAttempts) {
        std::vector<MovementResult> results;
        results.reserve(numAttempts);
        for (int i = 0; i < numAttempts; ++i) {
            results.push_back(attemptInsertion(state, params));
        }
        return results;
    }
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_INTERFACE_HPP
