#include "SimulationStatistics.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

SimulationStatistics::SimulationStatistics() 
    : lastUpdateTime_(std::chrono::steady_clock::now()) {
    // Initialize common move types
    moveStats_["insert"] = MoveStats();
    moveStats_["delete"] = MoveStats();
    moveStats_["translate"] = MoveStats();
    moveStats_["rotate"] = MoveStats();
}

void SimulationStatistics::recordMove(const std::string& moveType, 
                                     const std::string& fragmentName,
                                     bool accepted) {
    totalSteps_++;
    if (accepted) {
        totalAccepted_++;
    }
    
    // Update move-specific stats
    auto& stats = moveStats_[moveType];
    stats.attempts++;
    if (accepted) {
        stats.accepted++;
    }
    
    // Update fragment-specific stats
    auto& fragStats = fragmentStats_[fragmentName];
    fragStats.name = fragmentName;
    
    if (moveType == "insert") {
        fragStats.insertStats.attempts++;
        if (accepted) {
            fragStats.insertStats.accepted++;
            fragStats.currentCount++;
        }
    } else if (moveType == "delete") {
        fragStats.deleteStats.attempts++;
        if (accepted) {
            fragStats.deleteStats.accepted++;
            fragStats.currentCount--;
        }
    } else if (moveType == "translate") {
        fragStats.translateStats.attempts++;
        if (accepted) {
            fragStats.translateStats.accepted++;
        }
    } else if (moveType == "rotate") {
        fragStats.rotateStats.attempts++;
        if (accepted) {
            fragStats.rotateStats.accepted++;
        }
    }
}

void SimulationStatistics::recordEnergy(double energy) {
    currentEnergy_ = energy;
    energyHistory_.push_back(energy);
    energySum_ += energy;
    energySumSquared_ += energy * energy;
    energySamples_++;
}

void SimulationStatistics::recordFragmentCount(const std::string& name, int count) {
    fragmentStats_[name].currentCount = count;
}

void SimulationStatistics::recordStepTime(double seconds) {
    totalTime_ += seconds;
    stepTimes_.push_back(seconds);
}

void SimulationStatistics::updateStatistics() {
    // Update acceptance rates
    for (auto& [moveType, stats] : moveStats_) {
        stats.update();
    }
    
    for (auto& [fragName, fragStats] : fragmentStats_) {
        fragStats.insertStats.update();
        fragStats.deleteStats.update();
        fragStats.translateStats.update();
        fragStats.rotateStats.update();
    }
    
    updateAverages();
}

void SimulationStatistics::updateFragmentDensity(const std::string& name, double density) {
    fragmentStats_[name].density = density;
}

SimulationStatistics::MoveStats 
SimulationStatistics::getMoveStats(const std::string& moveType) const {
    auto it = moveStats_.find(moveType);
    return it != moveStats_.end() ? it->second : MoveStats();
}

SimulationStatistics::FragmentStats 
SimulationStatistics::getFragmentStats(const std::string& name) const {
    auto it = fragmentStats_.find(name);
    return it != fragmentStats_.end() ? it->second : FragmentStats();
}

double SimulationStatistics::getAverageEnergy() const {
    return energySamples_ > 0 ? energySum_ / energySamples_ : 0.0;
}

double SimulationStatistics::getEnergyStdDev() const {
    return calculateStdDev(energySum_, energySumSquared_, energySamples_);
}

double SimulationStatistics::getTotalAcceptanceRate() const {
    return totalSteps_ > 0 ? 
        static_cast<double>(totalAccepted_) / totalSteps_ : 0.0;
}

double SimulationStatistics::getStepsPerSecond() const {
    return totalTime_ > 0 ? totalSteps_ / totalTime_ : 0.0;
}

void SimulationStatistics::printSummary(int step) const {
    std::cout << "\n=== Step " << step << " Statistics ===" << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    
    // Overall statistics
    std::cout << "Total acceptance rate: " 
              << getTotalAcceptanceRate() * 100 << "%" << std::endl;
    std::cout << "Current energy: " << currentEnergy_ << " kJ/mol" << std::endl;
    std::cout << "Average energy: " << getAverageEnergy() 
              << " ± " << getEnergyStdDev() << " kJ/mol" << std::endl;
    
    // Move statistics
    std::cout << "\nMove acceptance rates:" << std::endl;
    for (const auto& [moveType, stats] : moveStats_) {
        if (stats.attempts > 0) {
            std::cout << "  " << moveType << ": " 
                     << stats.accepted << "/" << stats.attempts 
                     << " (" << stats.acceptanceRate * 100 << "%)" << std::endl;
        }
    }
    
    // Fragment statistics
    std::cout << "\nFragment counts:" << std::endl;
    for (const auto& [name, fragStats] : fragmentStats_) {
        std::cout << "  " << name << ": " << fragStats.currentCount;
        if (fragStats.density > 0) {
            std::cout << " (density: " << fragStats.density << " M)";
        }
        std::cout << std::endl;
    }
    
    // Performance
    std::cout << "\nPerformance: " << getStepsPerSecond() 
              << " steps/sec" << std::endl;
}

void SimulationStatistics::printDetailedStats() const {
    std::cout << "\n=== Detailed Statistics ===" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    
    for (const auto& [name, fragStats] : fragmentStats_) {
        std::cout << "\nFragment: " << name << std::endl;
        std::cout << "  Current count: " << fragStats.currentCount << std::endl;
        std::cout << "  Density: " << fragStats.density << " M" << std::endl;
        
        if (fragStats.insertStats.attempts > 0) {
            std::cout << "  Insert: " << fragStats.insertStats.accepted 
                     << "/" << fragStats.insertStats.attempts 
                     << " (" << fragStats.insertStats.acceptanceRate * 100 << "%)" 
                     << std::endl;
        }
        if (fragStats.deleteStats.attempts > 0) {
            std::cout << "  Delete: " << fragStats.deleteStats.accepted 
                     << "/" << fragStats.deleteStats.attempts 
                     << " (" << fragStats.deleteStats.acceptanceRate * 100 << "%)" 
                     << std::endl;
        }
        if (fragStats.translateStats.attempts > 0) {
            std::cout << "  Translate: " << fragStats.translateStats.accepted 
                     << "/" << fragStats.translateStats.attempts 
                     << " (" << fragStats.translateStats.acceptanceRate * 100 << "%)" 
                     << std::endl;
        }
        if (fragStats.rotateStats.attempts > 0) {
            std::cout << "  Rotate: " << fragStats.rotateStats.accepted 
                     << "/" << fragStats.rotateStats.attempts 
                     << " (" << fragStats.rotateStats.acceptanceRate * 100 << "%)" 
                     << std::endl;
        }
    }
}

std::string SimulationStatistics::formatStatistics() const {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2);
    
    ss << "Step " << totalSteps_ 
       << " | Acc: " << getTotalAcceptanceRate() * 100 << "%"
       << " | E: " << currentEnergy_ << " kJ/mol"
       << " | ";
    
    // Fragment counts
    for (const auto& [name, fragStats] : fragmentStats_) {
        ss << name << ": " << fragStats.currentCount << " ";
    }
    
    ss << "| " << getStepsPerSecond() << " steps/s";
    
    return ss.str();
}

void SimulationStatistics::reset() {
    totalSteps_ = 0;
    totalAccepted_ = 0;
    totalTime_ = 0.0;
    
    moveStats_.clear();
    fragmentStats_.clear();
    
    energyHistory_.clear();
    currentEnergy_ = 0.0;
    energySum_ = 0.0;
    energySumSquared_ = 0.0;
    energySamples_ = 0;
    
    stepTimes_.clear();
    lastUpdateTime_ = std::chrono::steady_clock::now();
}

void SimulationStatistics::updateAverages() {
    // Update timing
    auto now = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
        now - lastUpdateTime_);
    lastUpdateTime_ = now;
}

double SimulationStatistics::calculateStdDev(double sum, double sumSquared, int n) const {
    if (n <= 1) return 0.0;
    
    double mean = sum / n;
    double variance = (sumSquared / n) - (mean * mean);
    
    // Handle numerical errors
    if (variance < 0) variance = 0;
    
    return std::sqrt(variance);
}

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc