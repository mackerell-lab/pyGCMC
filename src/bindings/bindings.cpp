#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace pygcmc {
namespace bindings {

// Forward declarations of submodule initialization functions
void init_io(py::module& m);

PYBIND11_MODULE(pygcmc_python, m) {
    m.doc() = "Python bindings for PYGCMC library"; // Optional module docstring
    
    // Initialize submodules
    init_io(m);
}

} // namespace bindings
} // namespace pygcmc 