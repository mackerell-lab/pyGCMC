#pragma once

#include <atomic>
#include <stdexcept>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Global safety checks to prevent operations after cleanup
 */
class SafetyChecks {
private:
    static std::atomic<bool> cleanup_called;
    static std::atomic<int> initialization_count;
    
public:
    /**
     * @brief Mark that cleanup has been called
     */
    static void markCleanupCalled() {
        cleanup_called.store(true);
    }
    
    /**
     * @brief Mark that initialization has been done
     */
    static void markInitialized() {
        initialization_count.fetch_add(1);
        cleanup_called.store(false);
    }
    
    /**
     * @brief Check if it's safe to proceed with operations
     */
    static void checkSafeToOperate(const char* operation_name) {
        if (cleanup_called.load() && initialization_count.load() > 0) {
            throw std::runtime_error(
                std::string("Operation '") + operation_name + 
                "' called after cleanup without re-initialization. " +
                "Please re-initialize parameters before continuing."
            );
        }
    }
    
    /**
     * @brief Reset all safety flags (for testing only)
     */
    static void reset() {
        cleanup_called.store(false);
        initialization_count.store(0);
    }
};

// Initialize static members
inline std::atomic<bool> SafetyChecks::cleanup_called{false};
inline std::atomic<int> SafetyChecks::initialization_count{0};

/**
 * @brief RAII guard for operations
 */
class OperationGuard {
private:
    const char* operation_name;
    
public:
    explicit OperationGuard(const char* name) : operation_name(name) {
        SafetyChecks::checkSafeToOperate(operation_name);
    }
    
    ~OperationGuard() = default;
};

// Convenience macro for safety checks
#define SAFETY_CHECK(op) OperationGuard _guard(op)

} // namespace cpu
} // namespace platform
} // namespace pygcmc