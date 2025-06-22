#pragma once

#ifndef PYGCMC_MODEL_UTILS_TESTING_HPP
#define PYGCMC_MODEL_UTILS_TESTING_HPP

#include "UtilsFactory.hpp"
#include "UtilsValidation.hpp"
#include "../atom/AtomMain.hpp"
#include "../residue/ResidueMain.hpp"
#include "../molecule/MolecularMain.hpp"
#include <memory>
#include <string>
#include <iostream>
#include <chrono>
#include <vector>

namespace pygcmc {
namespace model {
namespace testing {

/**
 * @brief Run basic functionality tests for all components
 * @return true if all basic tests pass, false otherwise
 */
inline bool run_basic_tests() {
    try {
        // Test atom creation
        auto atom = std::make_shared<atom::Atom>(1, "CT1", "ALA", 1, "PROT", 0, 
                                                  0.0, 0.0, 0.0, 1.0, 12.01, 0.07, "C", false);
        if (!validation::validate_atom(*atom)) return false;
        
        // Test residue creation
        auto residue = factory::create_alanine_residue(1);
        if (!validation::validate_residue(*residue)) return false;
        
        // Test molecular system
        auto molecular = std::make_shared<molecule::Molecular>();
        molecular->add_residue(residue);
        if (!validation::validate_molecular_system(*molecular)) return false;
        
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

/**
 * @brief Test factory functions
 * @return true if all factory tests pass, false otherwise
 */
inline bool test_factory_functions() {
    try {
        // Test water molecule creation
        auto water = factory::create_water_molecule(1);
        if (!water || !validation::validate_residue(*water)) return false;
        
        // Test alanine residue creation
        auto ala = factory::create_alanine_residue(2);
        if (!ala || !validation::validate_residue(*ala)) return false;
        
        // Test molecular system validation
        auto molecular = std::make_shared<molecule::Molecular>();
        molecular->add_residue(water);
        molecular->add_residue(ala);
        
        if (!validation::validate_molecular_system(*molecular)) return false;
        
        // Test system summary
        std::string summary = validation::get_system_summary(*molecular);
        if (summary.empty()) return false;
        
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

/**
 * @brief Test simple system creation
 * @return true if simple system tests pass, false otherwise
 */
inline bool test_simple_system_creation() {
    try {
        // Test default simple system
        auto simple_system = factory::create_simple_system();
        if (!simple_system || !validation::validate_molecular_system(*simple_system)) {
            return false;
        }
        
        // Test custom simple system
        auto custom_system = factory::create_simple_system(5, 2);
        if (!custom_system || !validation::validate_molecular_system(*custom_system)) {
            return false;
        }
        
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

/**
 * @brief Test validation functions
 * @return true if all validation tests pass, false otherwise
 */
inline bool test_validation_functions() {
    try {
        // Create test system
        auto water = factory::create_water_molecule(1);
        auto ala = factory::create_alanine_residue(2);
        auto molecular = factory::create_molecular_system({water, ala});
        
        // Test individual validation functions
        if (!validation::validate_residue(*water)) return false;
        if (!validation::validate_residue(*ala)) return false;
        if (!validation::validate_molecular_system(*molecular)) return false;
        
        // Test detailed validation
        auto [is_valid, report] = validation::validate_system_detailed(*molecular);
        if (!is_valid || report.empty()) return false;
        
        // Test issue checking
        auto issues = validation::check_system_issues(*molecular);
        // Should have no serious issues for a well-formed system
        
        return true;
    } catch (const std::exception&) {
        return false;
    }
}

/**
 * @brief Run all tests
 * @return true if all tests pass, false otherwise
 */
inline bool run_all_tests() {
    return run_basic_tests() && 
           test_factory_functions() && 
           test_simple_system_creation() &&
           test_validation_functions();
}

/**
 * @brief Run tests with detailed output
 * @param verbose If true, print detailed information about each test
 * @return true if all tests pass, false otherwise
 */
inline bool run_tests_verbose(bool verbose = true) {
    bool all_passed = true;
    
    if (verbose) std::cout << "=== Running Model Utils Tests ===\n";
    
    // Test basic functionality
    if (verbose) std::cout << "Testing basic functionality... ";
    bool basic_ok = run_basic_tests();
    if (verbose) std::cout << (basic_ok ? "PASS" : "FAIL") << "\n";
    all_passed &= basic_ok;
    
    // Test factory functions
    if (verbose) std::cout << "Testing factory functions... ";
    bool factory_ok = test_factory_functions();
    if (verbose) std::cout << (factory_ok ? "PASS" : "FAIL") << "\n";
    all_passed &= factory_ok;
    
    // Test simple system creation
    if (verbose) std::cout << "Testing simple system creation... ";
    bool simple_ok = test_simple_system_creation();
    if (verbose) std::cout << (simple_ok ? "PASS" : "FAIL") << "\n";
    all_passed &= simple_ok;
    
    // Test validation functions
    if (verbose) std::cout << "Testing validation functions... ";
    bool validation_ok = test_validation_functions();
    if (verbose) std::cout << (validation_ok ? "PASS" : "FAIL") << "\n";
    all_passed &= validation_ok;
    
    if (verbose) {
        std::cout << "=== Test Results ===\n";
        std::cout << "Overall: " << (all_passed ? "PASS" : "FAIL") << "\n";
    }
    
    return all_passed;
}

/**
 * @brief Performance test for factory functions
 * @param num_iterations Number of iterations to run
 * @return Average time per iteration in milliseconds
 */
inline double performance_test_factory(int num_iterations = 1000) {
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < num_iterations; ++i) {
        auto water = factory::create_water_molecule(i);
        auto ala = factory::create_alanine_residue(i + 1000);
        auto system = factory::create_molecular_system({water, ala});
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    return static_cast<double>(duration.count()) / 1000.0 / num_iterations; // ms per iteration
}

} // namespace testing

/**
 * @brief Utilities namespace for backward compatibility
 */
namespace utils {

/**
 * @brief Run comprehensive module validation tests
 * @return true if all tests pass, false otherwise
 */
inline bool validate_module() {
    return testing::run_all_tests();
}

/**
 * @brief Run all factory function tests
 * @return true if factory tests pass, false otherwise
 */
inline bool test_all_factories() {
    return testing::test_factory_functions();
}

/**
 * @brief Run basic component tests
 * @return true if basic tests pass, false otherwise
 */
inline bool test_basic_components() {
    return testing::run_basic_tests();
}

/**
 * @brief Check if all module components are available
 * @return true if all components are working properly, false otherwise
 */
inline bool check_module_integrity() {
    try {
        bool factory_ok = testing::run_basic_tests();
        bool compatibility_ok = testing::test_factory_functions();
        
        return factory_ok && compatibility_ok;
    } catch (const std::exception&) {
        return false;
    }
}

} // namespace utils
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_UTILS_TESTING_HPP 