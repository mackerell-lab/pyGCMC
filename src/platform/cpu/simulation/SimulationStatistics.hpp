#ifndef PYGCMC_PLATFORM_CPU_SIMULATION_SIMULATION_STATISTICS_HPP
#define PYGCMC_PLATFORM_CPU_SIMULATION_SIMULATION_STATISTICS_HPP

#include <map>
#include <vector>
#include <string>
#include <chrono>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

/**
 * @brief Statistics tracking for GCMC simulation
 * 
 * This class manages all statistics collection and analysis
 * for the simulation, tracking moves, energy, and performance.
 */
class SimulationStatistics {
public:
    /**
     * @brief Move statistics structure
     */
    struct MoveStats {
        int attempts = 0;
        int accepted = 0;
        double acceptanceRate = 0.0;
        
        void update() {
            acceptanceRate = attempts > 0 ? 
                static_cast<double>(accepted) / attempts : 0.0;
        }
    };
    
    /**
     * @brief Fragment type statistics
     */
    struct FragmentStats {
        std::string name;
        int typeId;
        int currentCount = 0;
        double density = 0.0;
        MoveStats insertStats;
        MoveStats deleteStats;
        MoveStats translateStats;
        MoveStats rotateStats;
    };
    
    // Constructor
    SimulationStatistics();
    
    // Recording methods
    void recordMove(const std::string& moveType, 
                   const std::string& fragmentName,
                   bool accepted);
    void recordEnergy(double energy);
    void recordFragmentCount(const std::string& name, int count);
    void recordStepTime(double seconds);
    
    // Update methods
    void updateStatistics();
    void updateFragmentDensity(const std::string& name, double density);
    
    // Query methods
    MoveStats getMoveStats(const std::string& moveType) const;
    FragmentStats getFragmentStats(const std::string& name) const;
    double getAverageEnergy() const;
    double getEnergyStdDev() const;
    double getTotalAcceptanceRate() const;
    double getStepsPerSecond() const;
    
    // Output methods
    void printSummary(int step) const;
    void printDetailedStats() const;
    std::string formatStatistics() const;
    
    // Reset
    void reset();
    
    // Getters
    int getTotalSteps() const { return totalSteps_; }
    int getTotalAccepted() const { return totalAccepted_; }
    double getTotalTime() const { return totalTime_; }
    const std::vector<double>& getEnergyHistory() const { 
        return energyHistory_; 
    }
    
private:
    // Overall statistics
    int totalSteps_ = 0;
    int totalAccepted_ = 0;
    double totalTime_ = 0.0;  // seconds
    
    // Move-specific statistics
    std::map<std::string, MoveStats> moveStats_;
    
    // Fragment-specific statistics
    std::map<std::string, FragmentStats> fragmentStats_;
    
    // Energy tracking
    std::vector<double> energyHistory_;
    double currentEnergy_ = 0.0;
    double energySum_ = 0.0;
    double energySumSquared_ = 0.0;
    int energySamples_ = 0;
    
    // Timing
    std::chrono::steady_clock::time_point lastUpdateTime_;
    std::vector<double> stepTimes_;
    
    // Helper methods
    void updateAverages();
    double calculateStdDev(double sum, double sumSquared, int n) const;
};

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_SIMULATION_SIMULATION_STATISTICS_HPP