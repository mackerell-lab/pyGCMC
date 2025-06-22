#pragma once

#ifndef PYGCMC_MODEL_UTILS_MAIN_HPP
#define PYGCMC_MODEL_UTILS_MAIN_HPP

/**
 * @file UtilsMain.hpp
 * @brief Main utilities header for the Model module
 * 
 * This header provides a single point of access to all utility functions
 * for the Model module. It includes factory functions, validation utilities,
 * testing frameworks, and compatibility layers.
 * 
 * Usage:
 *   #include "utils/UtilsMain.hpp"
 *   
 * For specific functionality, you can include individual headers:
 *   #include "utils/UtilsFactory.hpp"      // Factory functions
 *   #include "utils/UtilsValidation.hpp"   // Validation utilities
 *   #include "utils/UtilsTesting.hpp"      // Testing framework
 *   #include "utils/UtilsInfo.hpp"         // Version and info
 *   #include "utils/UtilsCompatibility.hpp" // Backward compatibility
 */

// Core functionality headers
#include "UtilsFactory.hpp"
#include "UtilsValidation.hpp"
#include "UtilsTesting.hpp"
#include "UtilsInfo.hpp"
#include "UtilsCompatibility.hpp"

namespace pygcmc {
namespace model {

/**
 * @brief Main utilities namespace that aggregates all utility functions
 * 
 * This namespace provides convenient access to all utility functions
 * without needing to know the specific sub-namespace.
 */
namespace utils {

// === Factory Functions ===
using namespace factory;

// === Validation Functions ===
using namespace validation;

// === Testing Functions ===
using namespace testing;

// === Information Functions ===
using namespace info;

/**
 * @brief Quick utility initialization and self-test
 * @param verbose If true, print detailed test results
 * @return true if all utilities are working correctly
 */
inline bool initialize_and_test(bool verbose = false) {
    if (verbose) {
        std::cout << "=== Model Utils Initialization ===\n";
        std::cout << get_module_info() << "\n";
        std::cout << "Running self-tests...\n";
    }
    
    bool result = run_tests_verbose(verbose);
    
    if (verbose) {
        std::cout << "Initialization " << (result ? "successful" : "failed") << ".\n";
    }
    
    return result;
}

/**
 * @brief Quick system creation for testing and examples
 * @return A simple molecular system for demonstration
 */
inline auto create_demo_system() {
    return create_simple_system(5, 1); // 5 waters, 1 alanine
}

/**
 * @brief Comprehensive system analysis
 * @param system The molecular system to analyze
 * @return Detailed analysis report
 */
inline std::string analyze_system(const molecule::Molecular& system) {
    auto [is_valid, report] = validate_system_detailed(system);
    return report;
}

} // namespace utils

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_UTILS_MAIN_HPP 