#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_STATISTICS_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_STATISTICS_HPP

#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include "ProposalTypes.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

/**
 * P2: Enhanced proposal statistics for observability
 */
struct ProposalStatistics {
    // Basic counts per mode
    std::map<ProposalType, int> attempts;
    std::map<ProposalType, int> accepts;
    
    // Timing statistics (ms)
    std::vector<double> proposalTimes;
    std::vector<double> findCavTimes;
    
    // Mode distribution
    ProposalType currentMode = ProposalType::Uniform;
    int modeTransitions = 0;
    
    // Cavity statistics
    int lastNcav = 0;
    double lastOccupancy = 0.0;
    
    // Fallback counters
    int fallbackNoCavities = 0;       // Fallback due to no cavities found
    int fallbackTimeout = 0;          // Fallback due to timeout
    int fallbackInvalidMode = 0;      // Fallback due to invalid mode
    int autoSwitches = 0;             // Automatic mode switches
    std::map<std::string, int> switchReasons;  // Reasons for mode switches
    
    // Performance percentiles
    struct Percentiles {
        double p50 = 0.0;
        double p90 = 0.0;
        double p95 = 0.0;
        double p99 = 0.0;
    };
    
    // Helper to calculate percentiles
    Percentiles calculatePercentiles(const std::vector<double>& times) const {
        Percentiles p;
        if (times.empty()) return p;
        
        std::vector<double> sorted = times;
        std::sort(sorted.begin(), sorted.end());
        
        auto getPercentile = [&sorted](double percentile) {
            size_t idx = static_cast<size_t>(sorted.size() * percentile / 100.0);
            if (idx >= sorted.size()) idx = sorted.size() - 1;
            return sorted[idx];
        };
        
        p.p50 = getPercentile(50);
        p.p90 = getPercentile(90);
        p.p95 = getPercentile(95);
        p.p99 = getPercentile(99);
        
        return p;
    }
    
    // Get acceptance rate for a mode
    double getAcceptanceRate(ProposalType type) const {
        auto attemptIt = attempts.find(type);
        auto acceptIt = accepts.find(type);
        
        if (attemptIt == attempts.end() || attemptIt->second == 0) {
            return 0.0;
        }
        
        int acceptCount = (acceptIt != accepts.end()) ? acceptIt->second : 0;
        return static_cast<double>(acceptCount) / attemptIt->second;
    }
    
    // Get total attempts
    int getTotalAttempts() const {
        int total = 0;
        for (const auto& [type, count] : attempts) {
            total += count;
        }
        return total;
    }
    
    // Get total accepts
    int getTotalAccepts() const {
        int total = 0;
        for (const auto& [type, count] : accepts) {
            total += count;
        }
        return total;
    }
    
    // Record an attempt
    void recordAttempt(ProposalType type, bool accepted, 
                       double proposalTimeMs = -1.0, 
                       double findCavTimeMs = -1.0) {
        attempts[type]++;
        if (accepted) {
            accepts[type]++;
        }
        
        if (proposalTimeMs >= 0) {
            proposalTimes.push_back(proposalTimeMs);
            // Keep only last 1000 entries for memory efficiency
            if (proposalTimes.size() > 1000) {
                proposalTimes.erase(proposalTimes.begin());
            }
        }
        
        if (findCavTimeMs >= 0) {
            findCavTimes.push_back(findCavTimeMs);
            if (findCavTimes.size() > 1000) {
                findCavTimes.erase(findCavTimes.begin());
            }
        }
    }
    
    // Reset statistics
    void reset() {
        attempts.clear();
        accepts.clear();
        proposalTimes.clear();
        findCavTimes.clear();
        modeTransitions = 0;
        lastNcav = 0;
        lastOccupancy = 0.0;
        fallbackNoCavities = 0;
        fallbackTimeout = 0;
        fallbackInvalidMode = 0;
        autoSwitches = 0;
        switchReasons.clear();
    }
    
    // Record fallback event
    void recordFallback(const std::string& reason) {
        if (reason == "no_cavities") {
            fallbackNoCavities++;
        } else if (reason == "timeout") {
            fallbackTimeout++;
        } else if (reason == "invalid_mode") {
            fallbackInvalidMode++;
        }
        switchReasons[reason]++;
    }
    
    // Record mode switch
    void recordModeSwitch(ProposalType from, ProposalType to, const std::string& reason) {
        (void)from;  // Suppress unused parameter warning (could be used for logging)
        modeTransitions++;
        if (reason.find("auto") != std::string::npos) {
            autoSwitches++;
        }
        switchReasons[reason]++;
        currentMode = to;
    }
    
    // Get summary
    std::string getSummary() const {
        std::string summary = "ProposalStatistics:\n";
        summary += "  Total: " + std::to_string(getTotalAttempts()) + " attempts, ";
        summary += std::to_string(getTotalAccepts()) + " accepts\n";
        
        for (const auto& [type, count] : attempts) {
            double rate = getAcceptanceRate(type);
            summary += "  " + proposalTypeToString(type) + ": ";
            summary += std::to_string(count) + " attempts, ";
            summary += "accept rate = " + std::to_string(rate) + "\n";
        }
        
        if (!proposalTimes.empty()) {
            auto p = calculatePercentiles(proposalTimes);
            summary += "  Proposal time (ms): p50=" + std::to_string(p.p50);
            summary += ", p90=" + std::to_string(p.p90) + "\n";
        }
        
        if (!findCavTimes.empty()) {
            auto p = calculatePercentiles(findCavTimes);
            summary += "  FindCav time (ms): p50=" + std::to_string(p.p50);
            summary += ", p90=" + std::to_string(p.p90) + "\n";
        }
        
        // Fallback statistics
        if (fallbackNoCavities > 0 || fallbackTimeout > 0) {
            summary += "  Fallbacks: no_cavities=" + std::to_string(fallbackNoCavities);
            summary += ", timeout=" + std::to_string(fallbackTimeout) + "\n";
        }
        
        if (autoSwitches > 0) {
            summary += "  Auto switches: " + std::to_string(autoSwitches) + "\n";
        }
        
        if (!switchReasons.empty()) {
            summary += "  Switch reasons:\n";
            for (const auto& [reason, count] : switchReasons) {
                summary += "    " + reason + ": " + std::to_string(count) + "\n";
            }
        }
        
        return summary;
    }
    
private:
    std::string proposalTypeToString(ProposalType type) const {
        switch (type) {
            case ProposalType::Uniform: return "Uniform";
            case ProposalType::Cavity: return "Cavity";
            case ProposalType::Color: return "Color";
            case ProposalType::Cluster: return "Cluster";
            case ProposalType::Adaptive: return "Adaptive";
            default: return "Unknown";
        }
    }
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_PROPOSAL_STATISTICS_HPP