#pragma once

#ifndef PYGCMC_MODEL_UTILS_MODULE_UTILS_HPP
#define PYGCMC_MODEL_UTILS_MODULE_UTILS_HPP

#include "ModelFactory.hpp"

namespace pygcmc {
namespace model {
namespace utils {

/**
 * @brief Run comprehensive module validation tests
 */
inline bool validate_module() {
    return factory::run_all_tests();
}

/**
 * @brief Run all factory function tests
 */
inline bool test_all_factories() {
    return factory::test_factory_functions();
}

/**
 * @brief Run basic component tests
 */
inline bool test_basic_components() {
    return factory::run_basic_tests();
}

/**
 * @brief Check if all module components are available
 */
inline bool check_module_integrity() {
    try {
        // Test that all major components can be instantiated
        bool factory_ok = factory::run_basic_tests();
        bool compatibility_ok = factory::test_factory_functions();
        
        return factory_ok && compatibility_ok;
    } catch (const std::exception&) {
        return false;
    }
}

} // namespace utils
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_UTILS_MODULE_UTILS_HPP 