#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_STATISTICS_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_STATISTICS_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

/**
 * @brief Smart statistics collector for GCMC simulations
 * 
 * Automatically adjusts sampling frequency based on variance
 * to balance accuracy and performance
 */
class GCMCStatistics {
public:
    struct Stats {
        double mean = 0.0;
        double variance = 0.0;
        double min = 0.0;
        double max = 0.0;
        int count = 0;
    };
    
    struct Sample {
        int step;
        int particleCount;
        double energy;
        double acceptanceRate;
        double temperature;
        double activity;
    };
    
private:
    std::vector<Sample> samples_;
    int samplingInterval_ = 1000;  // Default sampling interval
    int minInterval_ = 10;         // Minimum interval (during equilibration)
    int maxInterval_ = 10000;      // Maximum interval (when equilibrated)
    bool autoAdjust_ = false;      // Auto-adjust sampling frequency
    int lastSampleStep_ = 0;
    
    // Running statistics
    double runningMean_ = 0.0;
    double runningM2_ = 0.0;  // For variance calculation
    int runningCount_ = 0;
    
    // Variance threshold for auto-adjustment
    double varianceThreshold_ = 0.01;
    
public:
    GCMCStatistics() = default;
    
    /**
     * Enable/disable auto-adjustment of sampling frequency
     */
    void setAutoAdjust(bool enable) {
        autoAdjust_ = enable;
    }
    
    /**
     * Set sampling interval (steps between samples)
     */
    void setSamplingInterval(int interval) {
        samplingInterval_ = std::max(1, interval);
    }
    
    /**
     * Get current sampling interval
     */
    int getSamplingInterval() const {
        return samplingInterval_;
    }
    
    /**
     * Set interval bounds for auto-adjustment
     */
    void setIntervalBounds(int minInt, int maxInt) {
        minInterval_ = std::max(1, minInt);
        maxInterval_ = std::max(minInterval_, maxInt);
    }
    
    /**
     * Check if we should sample at this step
     */
    bool shouldSample(int step) const {
        return (step - lastSampleStep_) >= samplingInterval_;
    }
    
    /**
     * Add a sample and optionally adjust sampling frequency
     */
    void addSample(int step, int particleCount, double energy, 
                   double acceptanceRate, double temperature, double activity) {
        // Record sample
        samples_.push_back({step, particleCount, energy, 
                          acceptanceRate, temperature, activity});
        lastSampleStep_ = step;
        
        // Update running statistics using Welford's algorithm
        runningCount_++;
        double delta = particleCount - runningMean_;
        runningMean_ += delta / runningCount_;
        double delta2 = particleCount - runningMean_;
        runningM2_ += delta * delta2;
        
        // Auto-adjust sampling frequency if enabled
        if (autoAdjust_ && runningCount_ > 10) {
            adjustSamplingFrequency();
        }
    }
    
    /**
     * Get statistics for particle count
     */
    Stats getParticleStats() const {
        Stats stats;
        if (samples_.empty()) return stats;
        
        std::vector<double> counts;
        for (const auto& s : samples_) {
            counts.push_back(s.particleCount);
        }
        
        stats.count = counts.size();
        stats.mean = std::accumulate(counts.begin(), counts.end(), 0.0) / stats.count;
        
        double sqSum = 0.0;
        for (double c : counts) {
            sqSum += (c - stats.mean) * (c - stats.mean);
        }
        stats.variance = sqSum / stats.count;
        
        auto minmax = std::minmax_element(counts.begin(), counts.end());
        stats.min = *minmax.first;
        stats.max = *minmax.second;
        
        return stats;
    }
    
    /**
     * Get statistics for energy
     */
    Stats getEnergyStats() const {
        Stats stats;
        if (samples_.empty()) return stats;
        
        std::vector<double> energies;
        for (const auto& s : samples_) {
            energies.push_back(s.energy);
        }
        
        stats.count = energies.size();
        stats.mean = std::accumulate(energies.begin(), energies.end(), 0.0) / stats.count;
        
        double sqSum = 0.0;
        for (double e : energies) {
            sqSum += (e - stats.mean) * (e - stats.mean);
        }
        stats.variance = sqSum / stats.count;
        
        auto minmax = std::minmax_element(energies.begin(), energies.end());
        stats.min = *minmax.first;
        stats.max = *minmax.second;
        
        return stats;
    }
    
    /**
     * Get all samples
     */
    const std::vector<Sample>& getSamples() const {
        return samples_;
    }
    
    /**
     * Clear all samples
     */
    void clear() {
        samples_.clear();
        runningMean_ = 0.0;
        runningM2_ = 0.0;
        runningCount_ = 0;
        lastSampleStep_ = 0;
    }
    
    /**
     * Get current variance (for auto-adjustment)
     */
    double getCurrentVariance() const {
        if (runningCount_ < 2) return 1.0;  // High variance initially
        return runningM2_ / (runningCount_ - 1);
    }
    
private:
    /**
     * Adjust sampling frequency based on variance
     */
    void adjustSamplingFrequency() {
        double variance = getCurrentVariance();
        double cv = std::sqrt(variance) / std::abs(runningMean_ + 1e-10);  // Coefficient of variation
        
        // If system is equilibrated (low CV), sample less frequently
        if (cv < 0.05) {
            samplingInterval_ = std::min(maxInterval_, 
                                        static_cast<int>(samplingInterval_ * 1.5));
        }
        // If system is changing rapidly (high CV), sample more frequently
        else if (cv > 0.20) {
            samplingInterval_ = std::max(minInterval_, 
                                        static_cast<int>(samplingInterval_ / 1.5));
        }
        // Otherwise keep current interval
    }
};

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_STATISTICS_HPP