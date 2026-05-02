// ============================================================================
// MovementStatsMethods.cpp - P2 enhanced statistics methods
// ============================================================================

#include "MovementMain.hpp"
#include "../bias/CavityBias.hpp"
#ifdef PYGCMC_USE_PROPOSAL_LAYER
#include "../proposal/ProposalMain.hpp"
#include "../proposal/ProposalStatistics.hpp"
#endif
#include <map>
#include <string>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

#ifdef PYGCMC_USE_PROPOSAL_LAYER
// Helper function to convert ProposalType to string
static std::string proposalTypeToString(ProposalType type) {
    switch (type) {
        case ProposalType::Uniform: return "uniform";
        case ProposalType::Cavity: return "cavity";
        case ProposalType::Color: return "color";
        case ProposalType::Cluster: return "cluster";
        case ProposalType::Adaptive: return "adaptive";
        default: return "unknown";
    }
}
#endif

/**
 * P2: Get proposal statistics as a map (for Python conversion)
 * When PYGCMC_USE_PROPOSAL_LAYER is enabled, returns real statistics.
 * Otherwise returns basic statistics from MovementModule.
 */
std::map<std::string, double> MovementModule::getProposalStatsMap() const {
    std::map<std::string, double> result;

#ifdef PYGCMC_USE_PROPOSAL_LAYER
    // If proposal layer is enabled and we have ProposalMain, get real stats
    if (proposalMain_) {
        const auto& stats = proposalMain_->getStatistics();

        // Basic statistics
        result["total_attempts"] = static_cast<double>(stats.getTotalAttempts());
        result["total_accepts"] = static_cast<double>(stats.getTotalAccepts());
        result["acceptance_rate"] = stats.getTotalAttempts() > 0 ?
            static_cast<double>(stats.getTotalAccepts()) / stats.getTotalAttempts() : 0.0;

        // Mode statistics
        result["current_mode"] = static_cast<double>(stats.currentMode);
        result["mode_transitions"] = static_cast<double>(stats.modeTransitions);
        result["auto_switches"] = static_cast<double>(stats.autoSwitches);

        // Per-mode statistics
        for (int i = 0; i <= 4; ++i) {
            auto type = static_cast<ProposalType>(i);
            std::string modeName = proposalTypeToString(type);
            result[modeName + "_attempts"] = static_cast<double>(
                stats.attempts.count(type) ? stats.attempts.at(type) : 0);
            result[modeName + "_accepts"] = static_cast<double>(
                stats.accepts.count(type) ? stats.accepts.at(type) : 0);
            result[modeName + "_accept_rate"] = stats.getAcceptanceRate(type);
        }

        // Timing percentiles
        if (!stats.proposalTimes.empty()) {
            auto p = stats.calculatePercentiles(stats.proposalTimes);
            result["proposal_time_p50_ms"] = p.p50;
            result["proposal_time_p90_ms"] = p.p90;
            result["proposal_time_p95_ms"] = p.p95;
            result["proposal_time_p99_ms"] = p.p99;
        } else {
            result["proposal_time_p50_ms"] = -1.0;
            result["proposal_time_p90_ms"] = -1.0;
            result["proposal_time_p95_ms"] = -1.0;
            result["proposal_time_p99_ms"] = -1.0;
        }

        if (!stats.findCavTimes.empty()) {
            auto p = stats.calculatePercentiles(stats.findCavTimes);
            result["findcav_time_p50_ms"] = p.p50;
            result["findcav_time_p90_ms"] = p.p90;
            result["findcav_time_p95_ms"] = p.p95;
            result["findcav_time_p99_ms"] = p.p99;
        } else {
            result["findcav_time_p50_ms"] = -1.0;
            result["findcav_time_p90_ms"] = -1.0;
            result["findcav_time_p95_ms"] = -1.0;
            result["findcav_time_p99_ms"] = -1.0;
        }

        // Fallback counters
        result["fallback_no_cavities"] = static_cast<double>(stats.fallbackNoCavities);
        result["fallback_timeout"] = static_cast<double>(stats.fallbackTimeout);
        result["fallback_invalid_mode"] = static_cast<double>(stats.fallbackInvalidMode);

        // Cavity info
        result["last_ncav"] = static_cast<double>(stats.lastNcav);
        result["last_occupancy"] = stats.lastOccupancy;
    } else {
        // ProposalMain not available, fall back to basic stats
        fillBasicProposalStats(result);
    }
#else
    // Proposal layer not enabled, return basic stats
    fillBasicProposalStats(result);
#endif

    // Always add cavity manager stats if available
    if (cavityManager_) {
        result["cavity_count"] = static_cast<double>(cavityManager_->getCavityCount());
        result["total_grid_points"] = static_cast<double>(cavityManager_->getTotalGridPoints());
    }

    return result;
}

// Helper to fill basic stats when proposal layer is not available
void MovementModule::fillBasicProposalStats(std::map<std::string, double>& result) const {
    // Basic statistics from MovementModule's own counters
    int totalAttempts = 0;
    int totalAccepts = 0;
    for (const auto& [moveType, stat] : stats_) {
        totalAttempts += stat.attempts;
        totalAccepts += stat.accepts;
    }

    result["total_attempts"] = static_cast<double>(totalAttempts);
    result["total_accepts"] = static_cast<double>(totalAccepts);
    result["acceptance_rate"] = totalAttempts > 0 ?
        static_cast<double>(totalAccepts) / totalAttempts : 0.0;

    // Current mode from params
    result["current_mode"] = static_cast<double>(params_.proposalMode);
    result["mode_transitions"] = 0.0;

    // Placeholders for detailed stats
    result["proposal_time_p50_ms"] = -1.0;
    result["proposal_time_p90_ms"] = -1.0;
    result["findcav_time_p50_ms"] = -1.0;
    result["findcav_time_p90_ms"] = -1.0;

    // Fallback counters (not tracked without proposal layer)
    result["fallback_no_cavities"] = 0.0;
    result["fallback_timeout"] = 0.0;
    result["fallback_invalid_mode"] = 0.0;
    result["auto_switches"] = 0.0;
}

/**
 * P2: Get cavity manager statistics as a map (for Python conversion)
 */
std::map<std::string, double> MovementModule::getCavityStatsMap() const {
    std::map<std::string, double> result;

    if (cavityManager_) {
        auto stats = cavityManager_->getStatistics();

        // Grid statistics
        result["total_grid_points"] = static_cast<double>(stats.totalGridPoints);
        result["occupied_points"] = static_cast<double>(stats.occupiedPoints);
        result["cavity_points"] = static_cast<double>(stats.cavityPoints);
        result["occupancy_ratio"] = stats.occupancyRatio;
        result["cavity_ratio"] = stats.cavityRatio;

        // Cache performance
        result["cache_hits"] = static_cast<double>(stats.cacheHits);
        result["cache_misses"] = static_cast<double>(stats.cacheMisses);
        double hit_rate = (stats.cacheHits + stats.cacheMisses) > 0 ?
            static_cast<double>(stats.cacheHits) / (stats.cacheHits + stats.cacheMisses) : 0.0;
        result["cache_hit_rate"] = hit_rate;

        // Cluster analysis
        result["cluster_count"] = static_cast<double>(stats.clusterCount);
        result["largest_cluster_size"] = static_cast<double>(stats.largestClusterSize);
        result["average_cluster_size"] = static_cast<double>(stats.averageClusterSize);

        // P2 enhanced metrics
        result["build_time_ms"] = stats.buildTimeMs;
        result["last_find_time_ms"] = stats.lastFindTimeMs;
        result["color_class_count"] = static_cast<double>(stats.colorClassCount);

        // Color class statistics
        result["color_min_cavities"] = static_cast<double>(stats.colorClassStats.minCavities);
        result["color_median_cavities"] = static_cast<double>(stats.colorClassStats.medianCavities);
        result["color_p95_cavities"] = static_cast<double>(stats.colorClassStats.p95Cavities);
        result["color_avg_cavities"] = stats.colorClassStats.avgCavities;

        // Incremental update stats
        result["incremental_updates"] = static_cast<double>(stats.incrementalUpdates);
        result["full_rebuilds"] = static_cast<double>(stats.fullRebuilds);
        result["dirty_ratio"] = stats.dirtyRatio;
    } else {
        // Return empty stats if no cavity manager
        result["total_grid_points"] = 0.0;
        result["cavity_points"] = 0.0;
        result["occupancy_ratio"] = 0.0;
    }

    return result;
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc
