#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_CONFIG_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_CONFIG_HPP

#include <cstdlib>
#include <string>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

/**
 * @brief Global configuration for GCMC simulations
 * 
 * Singleton pattern for managing runtime configuration
 * without impacting performance when features are disabled.
 */
class GCMCConfig {
public:
    static GCMCConfig& getInstance() {
        static GCMCConfig instance;
        return instance;
    }
    
    // Performance-related configuration
    struct Performance {
        bool storeProbability = false;      // Store acceptance probability
        bool enableDetailedStats = false;   // Enable detailed statistics
        int statsInterval = 1000;           // Statistics sampling interval
        bool useFastMath = true;           // Use fast math approximations
        bool enableCaching = true;         // Enable various caches
        int cacheSize = 1000;              // Max cache entries
    } performance;
    
    // Debug configuration
    struct Debug {
        bool verbose = false;               // Verbose output
        bool checkDetailedBalance = false;  // Check detailed balance
        bool trackMemory = false;          // Track memory usage
        bool logMoves = false;             // Log all moves
    } debug;
    
    // Advanced features
    struct Advanced {
        bool enableBatchOperations = false; // Batch insertions/deletions
        bool enableParallel = false;        // OpenMP parallelization
        int batchSize = 10;                // Default batch size
        int numThreads = 1;                // Number of OpenMP threads
    } advanced;
    
    // Update configuration from environment variables
    void updateFromEnvironment() {
        // Performance settings
        if (std::getenv("GCMC_STORE_PROB")) {
            performance.storeProbability = true;
        }
        if (std::getenv("GCMC_ENABLE_STATS")) {
            performance.enableDetailedStats = true;
            const char* interval = std::getenv("GCMC_STATS_INTERVAL");
            if (interval) {
                performance.statsInterval = std::atoi(interval);
            }
        }
        if (std::getenv("GCMC_NO_CACHE")) {
            performance.enableCaching = false;
        }
        
        // Debug settings
        if (std::getenv("GCMC_DEBUG")) {
            debug.verbose = true;
        }
        if (std::getenv("GCMC_CHECK_DB")) {
            debug.checkDetailedBalance = true;
        }
        if (std::getenv("GCMC_TRACK_MEMORY")) {
            debug.trackMemory = true;
        }
        
        // Advanced settings
        if (std::getenv("GCMC_BATCH_OPS")) {
            advanced.enableBatchOperations = true;
            const char* bsize = std::getenv("GCMC_BATCH_SIZE");
            if (bsize) {
                advanced.batchSize = std::atoi(bsize);
            }
        }
        
#ifdef PYGCMC_USE_OPENMP
        if (std::getenv("GCMC_PARALLEL")) {
            advanced.enableParallel = true;
            const char* nthreads = std::getenv("OMP_NUM_THREADS");
            if (nthreads) {
                advanced.numThreads = std::atoi(nthreads);
            }
        }
#endif
    }
    
    // Reset to defaults
    void reset() {
        performance = Performance();
        debug = Debug();
        advanced = Advanced();
    }
    
    // Check if any performance-impacting features are enabled
    bool hasPerformanceImpact() const {
        return performance.storeProbability || 
               performance.enableDetailedStats ||
               debug.verbose ||
               debug.checkDetailedBalance;
    }
    
private:
    GCMCConfig() {
        updateFromEnvironment();
    }
    
    // Prevent copying
    GCMCConfig(const GCMCConfig&) = delete;
    GCMCConfig& operator=(const GCMCConfig&) = delete;
};

// Convenience macros for cleaner code
#define GCMC_CONFIG GCMCConfig::getInstance()
#define GCMC_PERF_CONFIG GCMCConfig::getInstance().performance
#define GCMC_DEBUG_CONFIG GCMCConfig::getInstance().debug
#define GCMC_ADV_CONFIG GCMCConfig::getInstance().advanced

// Performance-critical inline checks - always check environment
inline bool shouldStoreProbability() {
    // Check environment variable each time to avoid state issues in parallel tests
    // This is still very fast as getenv is optimized
    return std::getenv("GCMC_STORE_PROB") != nullptr;
}

inline bool shouldCollectStats(int step) {
    // Check environment variable for consistency
    if (std::getenv("GCMC_ENABLE_STATS") == nullptr) {
        return false;
    }
    
    // Get interval from environment or use default
    const char* interval_str = std::getenv("GCMC_STATS_INTERVAL");
    int interval = interval_str ? std::atoi(interval_str) : 1000;
    
    return (step % interval == 0);
}

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_GCMC_CONFIG_HPP