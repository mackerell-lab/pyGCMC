// src/bindings/bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "BindingsModule.hpp"

namespace py = pybind11;

PYBIND11_MODULE(pygcmc, m) {
    m.doc() = "Python bindings for GCMC simulation library";
    
    // Initialize bindings - only io is refactored for now
    pygcmc::bindings::io::init_io_bindings(m);
    pygcmc::bindings::init_model(m);
    pygcmc::bindings::init_system(m);
    pygcmc::bindings::init_simulation_bindings(m);
} 