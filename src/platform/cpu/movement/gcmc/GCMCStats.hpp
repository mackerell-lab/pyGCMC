#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_STATS_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_STATS_HPP

#include <vector>
#include <map>
#include <string>
#include <chrono>
#include <fstream>
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

/**
 * @brief Comprehensive statistics collection for GCMC simulations
 *
 * This class tracks all relevant statistics during GCMC simulations:
 * - Move acceptance rates
 * - Energy distributions
 * - Molecule number fluctuations
 * - Chemical potential convergence
 * - Performance metrics
 */
class GCMCStats {
public:
    // Move statistics
    struct MoveStatistics {
        int attempts = 0;
        int accepted = 0;
        double totalEnergyChange = 0.0;
        double maxEnergyChange = 0.0;
        double minEnergyChange = 0.0;
        std::vector<double> energyHistory;

        double acceptanceRate() const {
            return attempts > 0 ? static_cast<double>(accepted) / attempts : 0.0;
        }

        double averageEnergyChange() const {
            return accepted > 0 ? totalEnergyChange / accepted : 0.0;
        }
    };

    // Fragment statistics
    struct FragmentStatistics {
        int typeId;
        std::string name;
        double chemicalPotential;

        // Number statistics
        std::vector<int> numberHistory;
        double averageNumber = 0.0;
        double variance = 0.0;
        int maxNumber = 0;
        int minNumber = 0;

        // Thermodynamic properties
        double effectiveChemicalPotential = 0.0;
        double fugacity = 0.0;
        double activity = 0.0;

        // Distribution analysis
        std::map<int, int> numberDistribution;

        double getFluctuation() const {
            return variance > 0 ? std::sqrt(variance) : 0.0;
        }
    };

    // Energy statistics
    struct EnergyStatistics {
        std::vector<double> totalEnergyHistory;
        std::vector<double> vdwEnergyHistory;
        std::vector<double> elecEnergyHistory;

        double averageTotal = 0.0;
        double averageVdw = 0.0;
        double averageElec = 0.0;

        double varianceTotal = 0.0;
        double varianceVdw = 0.0;
        double varianceElec = 0.0;

        // Energy distribution
        std::map<int, int> energyDistribution;
        double energyBinWidth = 1.0;  // kJ/mol

        double getHeatCapacity(double temperature) const {
            // Cv = (<E^2> - <E>^2) / (kT^2)
            double kT = 8.314e-3 * temperature;  // kJ/mol
            return varianceTotal / (kT * kT);
        }
    };

    // Constructor
    GCMCStats();
    ~GCMCStats();

    // Initialize
    void initialize(int nFragmentTypes);
    void setFragmentInfo(int typeId, const std::string& name, double mu);

    // Record move attempts
    void recordInsertion(int typeId, bool accepted, double deltaE);
    void recordDeletion(int typeId, bool accepted, double deltaE);
    void recordTranslation(int fragmentId, bool accepted, double deltaE);
    void recordRotation(int fragmentId, bool accepted, double deltaE);
    void recordSwap(int type1, int type2, bool accepted, double deltaE);

    // Record system state
    void recordState(const std::vector<int>& moleculeCounts,
                     double totalEnergy,
                     double vdwEnergy,
                     double elecEnergy);

    // Calculate averages
    void updateAverages();
    void calculateFluctuations();
    void calculateChemicalPotentials(double volume, double temperature);

    // Analysis
    double getCompressibility(double temperature, double volume) const;
    double getIsothermalCompressibility(int typeId, double temperature) const;
    std::pair<double, double> getEnergyFluctuation() const;

    // Convergence checks
    bool isNumberConverged(int typeId, double tolerance = 0.01) const;
    bool isEnergyConverged(double tolerance = 0.01) const;
    double getConvergenceMetric() const;

    // Output
    void print() const;
    void printDetailed() const;
    void saveToFile(const std::string& filename) const;
    void saveMoleculeDistribution(const std::string& filename, int typeId) const;
    void saveEnergyDistribution(const std::string& filename) const;

    // Time series analysis
    double getAutocorrelationTime(int typeId) const;
    double getEffectiveSampleSize(int typeId) const;
    std::vector<double> getBlockAverages(int typeId, int blockSize) const;

    // Access statistics
    const MoveStatistics& getInsertStats(int typeId) const;
    const MoveStatistics& getDeleteStats(int typeId) const;
    const MoveStatistics& getTranslateStats() const { return translateStats_; }
    const MoveStatistics& getRotateStats() const { return rotateStats_; }
    const FragmentStatistics& getFragmentStats(int typeId) const;
    const EnergyStatistics& getEnergyStats() const { return energyStats_; }

    // Performance metrics
    void recordStepTime(double timeMs);
    double getAverageStepTime() const { return averageStepTime_; }
    double getTotalTime() const;
    double getStepsPerSecond() const;

    // Reset
    void reset();
    void resetMoveStatistics();
    void clearHistory();

private:
    // Move statistics by type
    std::map<int, MoveStatistics> insertStats_;
    std::map<int, MoveStatistics> deleteStats_;
    MoveStatistics translateStats_;
    MoveStatistics rotateStats_;
    std::map<std::pair<int,int>, MoveStatistics> swapStats_;

    // Fragment statistics
    std::map<int, FragmentStatistics> fragmentStats_;

    // Energy statistics
    EnergyStatistics energyStats_;

    // Timing
    std::chrono::high_resolution_clock::time_point startTime_;
    double averageStepTime_;
    int totalSteps_;

    // Configuration
    int equilibrationSteps_;
    int currentStep_;

    // Helper methods
    double calculateAutocorrelation(const std::vector<double>& data, int lag) const;
    double calculateVariance(const std::vector<double>& data, double mean) const;
    void updateDistribution(std::map<int, int>& distribution,
                           double value, double binWidth);
};

/**
 * @brief Real-time statistics monitor
 */
class GCMCStatsMonitor {
public:
    GCMCStatsMonitor(GCMCStats* stats);

    // Start monitoring
    void startMonitoring(int intervalMs = 1000);
    void stopMonitoring();

    // Set alerts
    void setConvergenceAlert(double threshold);
    void setAcceptanceAlert(double minRate, double maxRate);

    // Check alerts
    bool checkAlerts();

private:
    GCMCStats* statistics_;
    bool monitoring_;
    int monitorInterval_;

    // Alert thresholds
    double convergenceThreshold_;
    double minAcceptance_;
    double maxAcceptance_;

    // Monitoring thread (if needed)
    void monitorLoop();
};

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_STATS_HPP
