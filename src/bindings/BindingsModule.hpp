#pragma once

/**
 * @file BindingsModule.hpp
 * @brief Unified entry point for all Python bindings
 * 
 * This module provides initialization functions for all Python binding functionality.
 * It follows the same pattern as SystemModule.hpp and IOModule.hpp for consistency.
 * 
 * **Architecture:**
 * - Modular architecture with clear separation by source module
 * - Each binding module is independently maintainable
 * - Well-structured: each file under 150 lines with focused responsibilities
 * - AI-friendly: optimal file sizes for AI processing and assistance
 * 
 * **Backward Compatibility:**
 * All original Python APIs are preserved. Existing Python code requires no changes.
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace pygcmc {
namespace bindings {

// Forward declarations for binding init functions
namespace io {
    void init_io_bindings(pybind11::module& m);
}

// Keep original function names for other modules (unchanged for now)
void init_model(pybind11::module& m);
void init_system(pybind11::module& m);  
void init_simulation_bindings(pybind11::module& m);

} // namespace bindings
} // namespace pygcmc

// Main pybind11 module definition
PYBIND11_MODULE(pygcmc, m) {
    m.doc() = "Python bindings for GCMC simulation library";
    
    // Initialize bindings - only io is refactored for now
    pygcmc::bindings::io::init_io_bindings(m);
    pygcmc::bindings::init_model(m);
    pygcmc::bindings::init_system(m);
    pygcmc::bindings::init_simulation_bindings(m);
}