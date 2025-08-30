#include "GCMCStats.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

// Constructor
GCMCStats::GCMCStats()
    : averageStepTime_(0.0),
      totalSteps_(0),
      equilibrationSteps_(0),
      currentStep_(0) {
    startTime_ = std::chrono::high_resolution_clock::now();
}

// Destructor
GCMCStats::~GCMCStats() {
}

// Initialize
void GCMCStats::initialize(int nFragmentTypes) {
    for (int i = 0; i < nFragmentTypes; ++i) {
        insertStats_[i] = MoveStatistics();
        deleteStats_[i] = MoveStatistics();
        fragmentStats_[i] = FragmentStatistics();
        fragmentStats_[i].typeId = i;
    }
    reset();
}

// Set fragment info
void GCMCStats::setFragmentInfo(int typeId, const std::string& name, double mu) {
    if (fragmentStats_.find(typeId) != fragmentStats_.end()) {
        fragmentStats_[typeId].name = name;
        fragmentStats_[typeId].chemicalPotential = mu;
    }
}

// Record insertion
void GCMCStats::recordInsertion(int typeId, bool accepted, double deltaE) {
    auto& stats = insertStats_[typeId];
    stats.attempts++;
    
    if (accepted) {
        stats.accepted++;
        stats.totalEnergyChange += deltaE;
        stats.maxEnergyChange = std::max(stats.maxEnergyChange, deltaE);
        stats.minEnergyChange = std::min(stats.minEnergyChange, deltaE);
        stats.energyHistory.push_back(deltaE);
    }
}

// Record deletion
void GCMCStats::recordDeletion(int typeId, bool accepted, double deltaE) {
    auto& stats = deleteStats_[typeId];
    stats.attempts++;
    
    if (accepted) {
        stats.accepted++;
        stats.totalEnergyChange += deltaE;
        stats.maxEnergyChange = std::max(stats.maxEnergyChange, deltaE);
        stats.minEnergyChange = std::min(stats.minEnergyChange, deltaE);
        stats.energyHistory.push_back(deltaE);
    }
}

// Record translation
void GCMCStats::recordTranslation(int fragmentId, bool accepted, double deltaE) {
    // Suppress unused parameter warning
    (void)fragmentId;
    
    translateStats_.attempts++;
    
    if (accepted) {
        translateStats_.accepted++;
        translateStats_.totalEnergyChange += deltaE;
        translateStats_.maxEnergyChange = std::max(translateStats_.maxEnergyChange, deltaE);
        translateStats_.minEnergyChange = std::min(translateStats_.minEnergyChange, deltaE);
        translateStats_.energyHistory.push_back(deltaE);
    }
}

// Record rotation
void GCMCStats::recordRotation(int fragmentId, bool accepted, double deltaE) {
    // Suppress unused parameter warning
    (void)fragmentId;
    
    rotateStats_.attempts++;
    
    if (accepted) {
        rotateStats_.accepted++;
        rotateStats_.totalEnergyChange += deltaE;
        rotateStats_.maxEnergyChange = std::max(rotateStats_.maxEnergyChange, deltaE);
        rotateStats_.minEnergyChange = std::min(rotateStats_.minEnergyChange, deltaE);
        rotateStats_.energyHistory.push_back(deltaE);
    }
}

// Record swap
void GCMCStats::recordSwap(int type1, int type2, bool accepted, double deltaE) {
    auto key = std::make_pair(type1, type2);
    auto& stats = swapStats_[key];
    stats.attempts++;
    
    if (accepted) {
        stats.accepted++;
        stats.totalEnergyChange += deltaE;
        stats.energyHistory.push_back(deltaE);
    }
}

// Record state
void GCMCStats::recordState(const std::vector<int>& moleculeCounts,
                           double totalEnergy,
                           double vdwEnergy,
                           double elecEnergy) {
    currentStep_++;
    
    // Update fragment counts
    for (size_t i = 0; i < moleculeCounts.size(); ++i) {
        if (fragmentStats_.find(i) != fragmentStats_.end()) {
            auto& fStats = fragmentStats_[i];
            fStats.numberHistory.push_back(moleculeCounts[i]);
            fStats.numberDistribution[moleculeCounts[i]]++;
            fStats.maxNumber = std::max(fStats.maxNumber, moleculeCounts[i]);
            fStats.minNumber = std::min(fStats.minNumber, moleculeCounts[i]);
        }
    }
    
    // Update energy
    energyStats_.totalEnergyHistory.push_back(totalEnergy);
    energyStats_.vdwEnergyHistory.push_back(vdwEnergy);
    energyStats_.elecEnergyHistory.push_back(elecEnergy);
    
    updateDistribution(energyStats_.energyDistribution, totalEnergy, 
                      energyStats_.energyBinWidth);
    
    totalSteps_++;
}

// Update averages
void GCMCStats::updateAverages() {
    // Fragment averages
    for (auto& [typeId, stats] : fragmentStats_) {
        if (!stats.numberHistory.empty()) {
            double sum = std::accumulate(stats.numberHistory.begin(), 
                                       stats.numberHistory.end(), 0.0);
            stats.averageNumber = sum / stats.numberHistory.size();
        }
    }
    
    // Energy averages
    if (!energyStats_.totalEnergyHistory.empty()) {
        double sumTotal = std::accumulate(energyStats_.totalEnergyHistory.begin(),
                                         energyStats_.totalEnergyHistory.end(), 0.0);
        energyStats_.averageTotal = sumTotal / energyStats_.totalEnergyHistory.size();
        
        double sumVdw = std::accumulate(energyStats_.vdwEnergyHistory.begin(),
                                       energyStats_.vdwEnergyHistory.end(), 0.0);
        energyStats_.averageVdw = sumVdw / energyStats_.vdwEnergyHistory.size();
        
        double sumElec = std::accumulate(energyStats_.elecEnergyHistory.begin(),
                                        energyStats_.elecEnergyHistory.end(), 0.0);
        energyStats_.averageElec = sumElec / energyStats_.elecEnergyHistory.size();
    }
}

// Calculate fluctuations
void GCMCStats::calculateFluctuations() {
    // Fragment number fluctuations
    for (auto& [typeId, stats] : fragmentStats_) {
        if (!stats.numberHistory.empty()) {
            stats.variance = calculateVariance(
                std::vector<double>(stats.numberHistory.begin(), stats.numberHistory.end()),
                stats.averageNumber);
        }
    }
    
    // Energy fluctuations
    if (!energyStats_.totalEnergyHistory.empty()) {
        energyStats_.varianceTotal = calculateVariance(energyStats_.totalEnergyHistory,
                                                       energyStats_.averageTotal);
        energyStats_.varianceVdw = calculateVariance(energyStats_.vdwEnergyHistory,
                                                     energyStats_.averageVdw);
        energyStats_.varianceElec = calculateVariance(energyStats_.elecEnergyHistory,
                                                      energyStats_.averageElec);
    }
}

// Calculate chemical potentials
void GCMCStats::calculateChemicalPotentials(double volume, double temperature) {
    double beta = 1.0 / (8.314e-3 * temperature);  // kJ/(mol*K)
    
    for (auto& [typeId, stats] : fragmentStats_) {
        if (stats.averageNumber > 0) {
            // μ_eff = kT * ln(<N> / V) - set chemical potential
            stats.effectiveChemicalPotential = 
                (1.0 / beta) * std::log(stats.averageNumber / volume);
            
            // Activity = <N> / V * exp(β * μ_set)
            stats.activity = (stats.averageNumber / volume) * 
                           std::exp(beta * stats.chemicalPotential);
            
            // Fugacity
            stats.fugacity = std::exp(beta * stats.effectiveChemicalPotential);
        }
    }
}

// Get compressibility
double GCMCStats::getCompressibility(double temperature, double volume) const {
    // κ_T = V * <(δN)^2> / (kT * <N>^2)
    double kT = 8.314e-3 * temperature;
    double totalVariance = 0.0;
    double totalNumber = 0.0;
    
    for (const auto& [typeId, stats] : fragmentStats_) {
        totalVariance += stats.variance;
        totalNumber += stats.averageNumber;
    }
    
    if (totalNumber > 0) {
        return volume * totalVariance / (kT * totalNumber * totalNumber);
    }
    
    return 0.0;
}

// Get isothermal compressibility
double GCMCStats::getIsothermalCompressibility(int typeId, double temperature) const {
    auto it = fragmentStats_.find(typeId);
    if (it == fragmentStats_.end()) return 0.0;
    
    const auto& stats = it->second;
    if (stats.averageNumber <= 0) return 0.0;
    
    double kT = 8.314e-3 * temperature;
    return stats.variance / (kT * stats.averageNumber);
}

// Get energy fluctuation
std::pair<double, double> GCMCStats::getEnergyFluctuation() const {
    return {energyStats_.averageTotal, std::sqrt(energyStats_.varianceTotal)};
}

// Check number convergence
bool GCMCStats::isNumberConverged(int typeId, double tolerance) const {
    auto it = fragmentStats_.find(typeId);
    if (it == fragmentStats_.end()) return false;
    
    const auto& stats = it->second;
    if (stats.numberHistory.size() < 100) return false;
    
    // Check last 10% vs previous 10%
    size_t n = stats.numberHistory.size();
    size_t blockSize = n / 10;
    
    double avgRecent = std::accumulate(stats.numberHistory.end() - blockSize,
                                      stats.numberHistory.end(), 0.0) / blockSize;
    double avgOld = std::accumulate(stats.numberHistory.begin(),
                                   stats.numberHistory.begin() + blockSize, 0.0) / blockSize;
    
    return std::abs(avgRecent - avgOld) / (avgOld + 1.0) < tolerance;
}

// Check energy convergence
bool GCMCStats::isEnergyConverged(double tolerance) const {
    if (energyStats_.totalEnergyHistory.size() < 100) return false;
    
    size_t n = energyStats_.totalEnergyHistory.size();
    size_t blockSize = n / 10;
    
    double avgRecent = std::accumulate(energyStats_.totalEnergyHistory.end() - blockSize,
                                      energyStats_.totalEnergyHistory.end(), 0.0) / blockSize;
    double avgOld = std::accumulate(energyStats_.totalEnergyHistory.begin(),
                                   energyStats_.totalEnergyHistory.begin() + blockSize, 0.0) / 
                                   blockSize;
    
    return std::abs(avgRecent - avgOld) / (std::abs(avgOld) + 1.0) < tolerance;
}

// Get convergence metric
double GCMCStats::getConvergenceMetric() const {
    double metric = 0.0;
    int count = 0;
    
    // Average convergence across all fragment types
    for (const auto& [typeId, stats] : fragmentStats_) {
        if (isNumberConverged(typeId, 0.01)) {
            metric += 1.0;
        }
        count++;
    }
    
    if (isEnergyConverged(0.01)) {
        metric += 1.0;
        count++;
    }
    
    return count > 0 ? metric / count : 0.0;
}

// Print statistics
void GCMCStats::print() const {
    std::cout << "\n=== GCMC Statistics ===" << std::endl;
    
    // Move statistics
    std::cout << "\n--- Move Acceptance Rates ---" << std::endl;
    for (const auto& [typeId, stats] : insertStats_) {
        std::cout << "Type " << typeId << " Insert: " 
                  << std::fixed << std::setprecision(3) 
                  << stats.acceptanceRate() << " (" 
                  << stats.accepted << "/" << stats.attempts << ")" << std::endl;
    }
    for (const auto& [typeId, stats] : deleteStats_) {
        std::cout << "Type " << typeId << " Delete: " 
                  << std::fixed << std::setprecision(3) 
                  << stats.acceptanceRate() << " (" 
                  << stats.accepted << "/" << stats.attempts << ")" << std::endl;
    }
    std::cout << "Translation: " << translateStats_.acceptanceRate() 
              << " (" << translateStats_.accepted << "/" 
              << translateStats_.attempts << ")" << std::endl;
    std::cout << "Rotation: " << rotateStats_.acceptanceRate() 
              << " (" << rotateStats_.accepted << "/" 
              << rotateStats_.attempts << ")" << std::endl;
    
    // Fragment statistics
    std::cout << "\n--- Fragment Statistics ---" << std::endl;
    for (const auto& [typeId, stats] : fragmentStats_) {
        std::cout << "Type " << typeId << " (" << stats.name << "):" << std::endl;
        std::cout << "  Average N: " << std::fixed << std::setprecision(2) 
                  << stats.averageNumber << " ± " << stats.getFluctuation() << std::endl;
        std::cout << "  Range: [" << stats.minNumber << ", " << stats.maxNumber << "]" 
                  << std::endl;
    }
    
    // Energy statistics
    std::cout << "\n--- Energy Statistics ---" << std::endl;
    std::cout << "Average Total: " << std::fixed << std::setprecision(2) 
              << energyStats_.averageTotal << " ± " 
              << std::sqrt(energyStats_.varianceTotal) << " kJ/mol" << std::endl;
    std::cout << "Average VdW: " << energyStats_.averageVdw << " kJ/mol" << std::endl;
    std::cout << "Average Elec: " << energyStats_.averageElec << " kJ/mol" << std::endl;
}

// Print detailed statistics
void GCMCStats::printDetailed() const {
    print();
    
    // Additional detailed information
    std::cout << "\n--- Thermodynamic Properties ---" << std::endl;
    for (const auto& [typeId, stats] : fragmentStats_) {
        std::cout << "Type " << typeId << ":" << std::endl;
        std::cout << "  Chemical Potential (set): " << stats.chemicalPotential 
                  << " kJ/mol" << std::endl;
        std::cout << "  Chemical Potential (eff): " << stats.effectiveChemicalPotential 
                  << " kJ/mol" << std::endl;
        std::cout << "  Activity: " << stats.activity << std::endl;
        std::cout << "  Fugacity: " << stats.fugacity << std::endl;
    }
    
    // Performance
    std::cout << "\n--- Performance ---" << std::endl;
    std::cout << "Total Steps: " << totalSteps_ << std::endl;
    std::cout << "Average Step Time: " << averageStepTime_ << " ms" << std::endl;
    std::cout << "Steps per Second: " << getStepsPerSecond() << std::endl;
    
    // Convergence
    std::cout << "\n--- Convergence ---" << std::endl;
    std::cout << "Overall Convergence: " << std::fixed << std::setprecision(1) 
              << getConvergenceMetric() * 100 << "%" << std::endl;
}

// Save to file
void GCMCStats::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    file << "# GCMC Statistics" << std::endl;
    file << "# Step\tN_total";
    for (const auto& [typeId, stats] : fragmentStats_) {
        file << "\tN_" << typeId;
    }
    file << "\tE_total\tE_vdw\tE_elec" << std::endl;
    
    // Write time series data
    size_t nSteps = energyStats_.totalEnergyHistory.size();
    for (size_t i = 0; i < nSteps; ++i) {
        file << i;
        
        int totalN = 0;
        for (const auto& [typeId, stats] : fragmentStats_) {
            if (i < stats.numberHistory.size()) {
                file << "\t" << stats.numberHistory[i];
                totalN += stats.numberHistory[i];
            } else {
                file << "\t0";
            }
        }
        
        file << "\t" << totalN;
        file << "\t" << energyStats_.totalEnergyHistory[i];
        
        if (i < energyStats_.vdwEnergyHistory.size()) {
            file << "\t" << energyStats_.vdwEnergyHistory[i];
        } else {
            file << "\t0";
        }
        
        if (i < energyStats_.elecEnergyHistory.size()) {
            file << "\t" << energyStats_.elecEnergyHistory[i];
        } else {
            file << "\t0";
        }
        
        file << std::endl;
    }
    
    file.close();
}

// Save molecule distribution
void GCMCStats::saveMoleculeDistribution(const std::string& filename, int typeId) const {
    auto it = fragmentStats_.find(typeId);
    if (it == fragmentStats_.end()) {
        std::cerr << "Error: Type " << typeId << " not found" << std::endl;
        return;
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    const auto& stats = it->second;
    
    file << "# Molecule Number Distribution for Type " << typeId 
         << " (" << stats.name << ")" << std::endl;
    file << "# N\tCount\tProbability" << std::endl;
    
    int totalCount = 0;
    for (const auto& [n, count] : stats.numberDistribution) {
        totalCount += count;
    }
    
    for (const auto& [n, count] : stats.numberDistribution) {
        file << n << "\t" << count << "\t" 
             << static_cast<double>(count) / totalCount << std::endl;
    }
    
    file.close();
}

// Save energy distribution
void GCMCStats::saveEnergyDistribution(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    
    file << "# Energy Distribution" << std::endl;
    file << "# E(kJ/mol)\tCount\tProbability" << std::endl;
    
    int totalCount = 0;
    for (const auto& [bin, count] : energyStats_.energyDistribution) {
        totalCount += count;
    }
    
    for (const auto& [bin, count] : energyStats_.energyDistribution) {
        double energy = bin * energyStats_.energyBinWidth;
        file << energy << "\t" << count << "\t" 
             << static_cast<double>(count) / totalCount << std::endl;
    }
    
    file.close();
}

// Get autocorrelation time
double GCMCStats::getAutocorrelationTime(int typeId) const {
    auto it = fragmentStats_.find(typeId);
    if (it == fragmentStats_.end()) return 0.0;
    
    const auto& history = it->second.numberHistory;
    if (history.size() < 10) return 0.0;
    
    std::vector<double> data(history.begin(), history.end());
    
    double tau = 0.0;
    for (int lag = 1; lag < static_cast<int>(data.size()) / 2; ++lag) {
        double corr = calculateAutocorrelation(data, lag);
        if (corr < 0.1) break;  // Cutoff when correlation becomes small
        tau += corr;
    }
    
    return 1.0 + 2.0 * tau;
}

// Get effective sample size
double GCMCStats::getEffectiveSampleSize(int typeId) const {
    double tau = getAutocorrelationTime(typeId);
    if (tau <= 0) return 0.0;
    
    auto it = fragmentStats_.find(typeId);
    if (it == fragmentStats_.end()) return 0.0;
    
    return it->second.numberHistory.size() / tau;
}

// Get block averages
std::vector<double> GCMCStats::getBlockAverages(int typeId, int blockSize) const {
    std::vector<double> blockAvgs;
    
    auto it = fragmentStats_.find(typeId);
    if (it == fragmentStats_.end()) return blockAvgs;
    
    const auto& history = it->second.numberHistory;
    
    for (size_t i = 0; i + blockSize <= history.size(); i += blockSize) {
        double sum = std::accumulate(history.begin() + i, 
                                    history.begin() + i + blockSize, 0.0);
        blockAvgs.push_back(sum / blockSize);
    }
    
    return blockAvgs;
}

// Get insert stats
const GCMCStats::MoveStatistics& GCMCStats::getInsertStats(int typeId) const {
    static const MoveStatistics empty;
    auto it = insertStats_.find(typeId);
    return (it != insertStats_.end()) ? it->second : empty;
}

// Get delete stats
const GCMCStats::MoveStatistics& GCMCStats::getDeleteStats(int typeId) const {
    static const MoveStatistics empty;
    auto it = deleteStats_.find(typeId);
    return (it != deleteStats_.end()) ? it->second : empty;
}

// Get fragment stats
const GCMCStats::FragmentStatistics& GCMCStats::getFragmentStats(int typeId) const {
    static const FragmentStatistics empty{};
    auto it = fragmentStats_.find(typeId);
    return (it != fragmentStats_.end()) ? it->second : empty;
}

// Record step time
void GCMCStats::recordStepTime(double timeMs) {
    averageStepTime_ = (averageStepTime_ * (totalSteps_ - 1) + timeMs) / totalSteps_;
}

// Get total time
double GCMCStats::getTotalTime() const {
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = now - startTime_;
    return elapsed.count();
}

// Get steps per second
double GCMCStats::getStepsPerSecond() const {
    double totalTimeS = getTotalTime() / 1000.0;
    return totalTimeS > 0 ? totalSteps_ / totalTimeS : 0.0;
}

// Reset
void GCMCStats::reset() {
    resetMoveStatistics();
    clearHistory();
    
    totalSteps_ = 0;
    currentStep_ = 0;
    averageStepTime_ = 0.0;
    startTime_ = std::chrono::high_resolution_clock::now();
}

// Reset move statistics
void GCMCStats::resetMoveStatistics() {
    for (auto& [typeId, stats] : insertStats_) {
        stats = MoveStatistics();
    }
    for (auto& [typeId, stats] : deleteStats_) {
        stats = MoveStatistics();
    }
    translateStats_ = MoveStatistics();
    rotateStats_ = MoveStatistics();
    swapStats_.clear();
}

// Clear history
void GCMCStats::clearHistory() {
    for (auto& [typeId, stats] : fragmentStats_) {
        stats.numberHistory.clear();
        stats.numberDistribution.clear();
        stats.averageNumber = 0.0;
        stats.variance = 0.0;
    }
    
    energyStats_.totalEnergyHistory.clear();
    energyStats_.vdwEnergyHistory.clear();
    energyStats_.elecEnergyHistory.clear();
    energyStats_.energyDistribution.clear();
    
    energyStats_.averageTotal = 0.0;
    energyStats_.averageVdw = 0.0;
    energyStats_.averageElec = 0.0;
    energyStats_.varianceTotal = 0.0;
    energyStats_.varianceVdw = 0.0;
    energyStats_.varianceElec = 0.0;
}

// Calculate autocorrelation
double GCMCStats::calculateAutocorrelation(const std::vector<double>& data, int lag) const {
    if (lag >= static_cast<int>(data.size())) return 0.0;
    
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    
    double numerator = 0.0;
    double denominator = 0.0;
    
    for (size_t i = 0; i < data.size() - lag; ++i) {
        numerator += (data[i] - mean) * (data[i + lag] - mean);
    }
    
    for (size_t i = 0; i < data.size(); ++i) {
        denominator += (data[i] - mean) * (data[i] - mean);
    }
    
    return denominator > 0 ? numerator / denominator : 0.0;
}

// Calculate variance
double GCMCStats::calculateVariance(const std::vector<double>& data, double mean) const {
    if (data.empty()) return 0.0;
    
    double sum = 0.0;
    for (double value : data) {
        double diff = value - mean;
        sum += diff * diff;
    }
    
    return sum / data.size();
}

// Update distribution
void GCMCStats::updateDistribution(std::map<int, int>& distribution,
                                  double value, double binWidth) {
    int bin = static_cast<int>(value / binWidth);
    distribution[bin]++;
}

// ============================================================================
// GCMCStatsMonitor Implementation
// ============================================================================

GCMCStatsMonitor::GCMCStatsMonitor(GCMCStats* stats)
    : statistics_(stats),
      monitoring_(false),
      monitorInterval_(1000),
      convergenceThreshold_(0.9),
      minAcceptance_(0.1),
      maxAcceptance_(0.5) {
}

// Start monitoring
void GCMCStatsMonitor::startMonitoring(int intervalMs) {
    monitoring_ = true;
    monitorInterval_ = intervalMs;
    // In practice, would start a monitoring thread
}

// Stop monitoring
void GCMCStatsMonitor::stopMonitoring() {
    monitoring_ = false;
}

// Set convergence alert
void GCMCStatsMonitor::setConvergenceAlert(double threshold) {
    convergenceThreshold_ = threshold;
}

// Set acceptance alert
void GCMCStatsMonitor::setAcceptanceAlert(double minRate, double maxRate) {
    minAcceptance_ = minRate;
    maxAcceptance_ = maxRate;
}

// Check alerts
bool GCMCStatsMonitor::checkAlerts() {
    if (!statistics_) return false;
    
    // Check convergence
    if (statistics_->getConvergenceMetric() < convergenceThreshold_) {
        std::cout << "Alert: Convergence below threshold" << std::endl;
        return true;
    }
    
    // Check acceptance rates
    const auto& translateStats = statistics_->getTranslateStats();
    if (translateStats.acceptanceRate() < minAcceptance_ ||
        translateStats.acceptanceRate() > maxAcceptance_) {
        std::cout << "Alert: Translation acceptance rate out of range" << std::endl;
        return true;
    }
    
    return false;
}

// Monitor loop
void GCMCStatsMonitor::monitorLoop() {
    while (monitoring_) {
        // Check and report statistics
        checkAlerts();
        
        // Sleep for interval
        // std::this_thread::sleep_for(std::chrono::milliseconds(monitorInterval_));
    }
}

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc