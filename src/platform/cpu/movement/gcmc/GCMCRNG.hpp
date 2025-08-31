#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_RNG_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_RNG_HPP

#include <random>
#include <memory>
#include <vector>
#include <cstdint>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

/**
 * @brief Centralized Random Number Generator for GCMC
 * 
 * This class provides a single source of randomness for all GCMC components,
 * ensuring reproducibility and avoiding multiple independent RNG instances.
 * 
 * Features:
 * - Single master RNG with controlled seeding
 * - Sub-stream creation for parallel components
 * - Reproducibility support with state recording
 * - Thread-safe access patterns (when used correctly)
 */
class GCMCRNG {
public:
    // Reproducibility record for debugging
    struct ReproRecord {
        uint64_t step;
        uint64_t stateHash;
        std::string component;
        
        ReproRecord(uint64_t s, uint64_t h, const std::string& c)
            : step(s), stateHash(h), component(c) {}
    };
    
    /**
     * @brief Constructor with explicit seed
     * @param seed Master seed for reproducibility
     * @param recordReproducibility If true, record state for replay
     */
    explicit GCMCRNG(uint64_t seed = 12345, bool recordReproducibility = false)
        : masterRng_(seed),
          currentStep_(0),
          recordRepro_(recordReproducibility),
          masterSeed_(seed) {
        
        // Initialize distributions
        uniform_ = std::uniform_real_distribution<double>(0.0, 1.0);
        normal_ = std::normal_distribution<double>(0.0, 1.0);
    }
    
    /**
     * @brief Get the master RNG engine
     * @warning Direct access - use with caution
     */
    std::mt19937_64& engine() { return masterRng_; }
    
    /**
     * @brief Generate uniform random number [0, 1)
     */
    double uniform() {
        recordAccess("uniform");
        return uniform_(masterRng_);
    }
    
    /**
     * @brief Generate uniform random number in range [min, max)
     */
    double uniform(double min, double max) {
        recordAccess("uniform_range");
        return min + (max - min) * uniform_(masterRng_);
    }
    
    /**
     * @brief Generate uniform integer in range [min, max]
     */
    int uniformInt(int min, int max) {
        recordAccess("uniform_int");
        std::uniform_int_distribution<int> dist(min, max);
        return dist(masterRng_);
    }
    
    /**
     * @brief Generate normal distributed random number
     */
    double normal(double mean = 0.0, double stddev = 1.0) {
        recordAccess("normal");
        return mean + stddev * normal_(masterRng_);
    }
    
    /**
     * @brief Create a sub-stream RNG for a component
     * Uses jump-ahead or split-mix seeding to avoid correlation
     */
    std::unique_ptr<std::mt19937_64> createSubStream(const std::string& componentName) {
        // Use splitmix64 to generate uncorrelated seed
        uint64_t subSeed = splitmix64(masterSeed_ ^ std::hash<std::string>{}(componentName));
        recordAccess("create_substream:" + componentName);
        return std::make_unique<std::mt19937_64>(subSeed);
    }
    
    /**
     * @brief Reset to initial state with same seed
     */
    void reset() {
        masterRng_.seed(masterSeed_);
        currentStep_ = 0;
        reproRecords_.clear();
    }
    
    /**
     * @brief Reset with new seed
     */
    void reseed(uint64_t newSeed) {
        masterSeed_ = newSeed;
        masterRng_.seed(newSeed);
        currentStep_ = 0;
        reproRecords_.clear();
    }
    
    /**
     * @brief Advance step counter (for tracking)
     */
    void advanceStep() { currentStep_++; }
    
    /**
     * @brief Get current step
     */
    uint64_t getCurrentStep() const { return currentStep_; }
    
    /**
     * @brief Get master seed
     */
    uint64_t getMasterSeed() const { return masterSeed_; }
    
    /**
     * @brief Get reproducibility records
     */
    const std::vector<ReproRecord>& getReproRecords() const { return reproRecords_; }
    
    /**
     * @brief Save RNG state for checkpoint
     */
    std::vector<uint64_t> saveState() const {
        std::vector<uint64_t> state;
        state.push_back(masterSeed_);
        state.push_back(currentStep_);
        // Note: std::mt19937_64 state is complex; for full state save,
        // would need to serialize internal state array
        return state;
    }
    
    /**
     * @brief Restore RNG state from checkpoint
     */
    void restoreState(const std::vector<uint64_t>& state) {
        if (state.size() >= 2) {
            masterSeed_ = state[0];
            currentStep_ = state[1];
            // Re-seed with combined seed to get deterministic state
            // Use splitmix64 to combine seed with step count
            uint64_t combinedSeed = splitmix64(masterSeed_ + currentStep_);
            masterRng_.seed(combinedSeed);
        }
    }
    
    /**
     * @brief Enable/disable reproducibility recording
     */
    void setRecordReproducibility(bool record) { recordRepro_ = record; }
    
private:
    // Master RNG engine
    std::mt19937_64 masterRng_;
    
    // Distributions
    std::uniform_real_distribution<double> uniform_;
    std::normal_distribution<double> normal_;
    
    // State tracking
    uint64_t masterSeed_;
    uint64_t currentStep_;
    
    // Reproducibility tracking
    bool recordRepro_;
    std::vector<ReproRecord> reproRecords_;
    
    /**
     * @brief Record access for reproducibility debugging
     */
    void recordAccess(const std::string& operation) {
        if (recordRepro_) {
            // Use current step and seed for hash - don't consume random numbers
            uint64_t stateHash = masterSeed_ ^ currentStep_ ^ std::hash<std::string>{}(operation);
            reproRecords_.emplace_back(currentStep_, stateHash, operation);
        }
    }
    
    /**
     * @brief SplitMix64 for generating uncorrelated seeds
     */
    static uint64_t splitmix64(uint64_t x) {
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27;
        x *= 0x94d049bb133111ebULL;
        x ^= x >> 31;
        return x;
    }
};

/**
 * @brief Global GCMC RNG singleton (optional pattern)
 * 
 * Usage:
 *   auto& rng = GlobalGCMCRNG::getInstance();
 *   double r = rng.uniform();
 */
class GlobalGCMCRNG {
public:
    static GCMCRNG& getInstance() {
        static GCMCRNG instance;
        return instance;
    }
    
    // Delete copy/move constructors
    GlobalGCMCRNG(const GlobalGCMCRNG&) = delete;
    GlobalGCMCRNG& operator=(const GlobalGCMCRNG&) = delete;
    
private:
    GlobalGCMCRNG() = default;
};

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_RNG_HPP