#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_STATISTICS_WRAPPER_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_STATISTICS_WRAPPER_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../proposal/ProposalStatistics.hpp"
#include "../bias/CavityBias.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

/**
 * P2: Wrapper to expose statistics to Python
 */
class MovementStatisticsWrapper {
public:
    /**
     * Convert ProposalStatistics to Python dict
     */
    static py::dict proposalStatsToDict(const ProposalStatistics& stats) {
        py::dict result;
        
        // Basic counts
        result["total_attempts"] = stats.getTotalAttempts();
        result["total_accepts"] = stats.getTotalAccepts();
        result["acceptance_rate"] = stats.getTotalAccepts() > 0 ? 
            static_cast<double>(stats.getTotalAccepts()) / stats.getTotalAttempts() : 0.0;
        
        // Per-mode statistics
        py::dict mode_stats;
        for (int i = 0; i <= 4; ++i) {
            ProposalType type = static_cast<ProposalType>(i);
            py::dict mode_info;
            
            auto attempts_it = stats.attempts.find(type);
            auto accepts_it = stats.accepts.find(type);
            
            int attempts = (attempts_it != stats.attempts.end()) ? attempts_it->second : 0;
            int accepts = (accepts_it != stats.accepts.end()) ? accepts_it->second : 0;
            
            mode_info["attempts"] = attempts;
            mode_info["accepts"] = accepts;
            mode_info["accept_rate"] = stats.getAcceptanceRate(type);
            
            mode_stats[proposalTypeToString(type)] = mode_info;
        }
        result["modes"] = mode_stats;
        
        // Current state
        result["current_mode"] = proposalTypeToString(stats.currentMode);
        result["mode_transitions"] = stats.modeTransitions;
        result["last_ncav"] = stats.lastNcav;
        result["last_occupancy"] = stats.lastOccupancy;
        
        // Timing percentiles
        if (!stats.proposalTimes.empty()) {
            auto p = stats.calculatePercentiles(stats.proposalTimes);
            py::dict timing;
            timing["p50"] = p.p50;
            timing["p90"] = p.p90;
            timing["p95"] = p.p95;
            timing["p99"] = p.p99;
            result["proposal_time_ms"] = timing;
        }
        
        if (!stats.findCavTimes.empty()) {
            auto p = stats.calculatePercentiles(stats.findCavTimes);
            py::dict timing;
            timing["p50"] = p.p50;
            timing["p90"] = p.p90;
            timing["p95"] = p.p95;
            timing["p99"] = p.p99;
            result["findcav_time_ms"] = timing;
        }
        
        return result;
    }
    
    /**
     * Convert CavityManager::Statistics to Python dict
     */
    static py::dict cavityStatsToDict(const CavityManager::Statistics& stats) {
        py::dict result;
        
        // Grid statistics
        result["total_grid_points"] = stats.totalGridPoints;
        result["occupied_points"] = stats.occupiedPoints;
        result["cavity_points"] = stats.cavityPoints;
        result["occupancy_ratio"] = stats.occupancyRatio;
        result["cavity_ratio"] = stats.cavityRatio;
        
        // Cache performance
        result["cache_hits"] = stats.cacheHits;
        result["cache_misses"] = stats.cacheMisses;
        double hit_rate = (stats.cacheHits + stats.cacheMisses) > 0 ?
            static_cast<double>(stats.cacheHits) / (stats.cacheHits + stats.cacheMisses) : 0.0;
        result["cache_hit_rate"] = hit_rate;
        
        // Cluster analysis
        result["cluster_count"] = stats.clusterCount;
        result["largest_cluster_size"] = stats.largestClusterSize;
        result["average_cluster_size"] = stats.averageClusterSize;
        
        // P2: Enhanced metrics
        result["build_time_ms"] = stats.buildTimeMs;
        result["last_find_time_ms"] = stats.lastFindTimeMs;
        result["color_class_count"] = stats.colorClassCount;
        
        // Color class statistics
        py::dict color_stats;
        color_stats["min_cavities"] = stats.colorClassStats.minCavities;
        color_stats["median_cavities"] = stats.colorClassStats.medianCavities;
        color_stats["p95_cavities"] = stats.colorClassStats.p95Cavities;
        color_stats["avg_cavities"] = stats.colorClassStats.avgCavities;
        result["color_class_stats"] = color_stats;
        
        // Incremental update stats (P3 prep)
        result["incremental_updates"] = stats.incrementalUpdates;
        result["full_rebuilds"] = stats.fullRebuilds;
        result["dirty_ratio"] = stats.dirtyRatio;
        
        return result;
    }
    
private:
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
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_STATISTICS_WRAPPER_HPP