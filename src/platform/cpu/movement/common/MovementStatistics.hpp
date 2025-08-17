#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_STATISTICS_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_STATISTICS_HPP

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

/**
 * Statistics tracking for movement operations
 */
struct MovementStatistics {
    int attempts = 0;
    int accepts = 0;
    double totalEnergyChange = 0.0;
    
    double acceptanceRate() const {
        return attempts > 0 ? static_cast<double>(accepts) / attempts : 0.0;
    }
    
    void reset() {
        attempts = 0;
        accepts = 0;
        totalEnergyChange = 0.0;
    }
    
    void recordAttempt(bool accepted, double energyChange = 0.0) {
        attempts++;
        if (accepted) {
            accepts++;
            totalEnergyChange += energyChange;
        }
    }
};

// Alias for compatibility
using Statistics = MovementStatistics;

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_STATISTICS_HPP