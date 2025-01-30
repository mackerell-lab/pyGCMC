// src/bindings/bindings.cpp

#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace pygcmc {
namespace bindings {

// Forward declarations of submodule initialization functions
void init_model(py::module& m);
void init_io(py::module& m);

PYBIND11_MODULE(pygcmc, m) {
    m.doc() = "Python bindings for PYGCMC library"; // Optional module docstring
    
    // Initialize submodules
    init_model(m);
    init_io(m);
}

} // namespace bindings
} // namespace pygcmc 