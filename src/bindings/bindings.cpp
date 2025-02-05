// src/bindings/bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace pygcmc {
namespace bindings {

// Forward declarations of submodule initialization functions
void init_model(py::module& m);
void init_io(py::module& m);
void init_system(py::module& m);
void init_simulation_bindings(py::module& m);

PYBIND11_MODULE(pygcmc, m) {
    m.doc() = "Python bindings for GCMC simulation library"; // optional module docstring
    
    // Initialize submodules
    init_model(m);
    init_io(m);
    init_system(m);
    init_simulation_bindings(m);
}

} // namespace bindings
} // namespace pygcmc 