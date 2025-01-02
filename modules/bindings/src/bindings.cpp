// modules/bindings/src/bindings.cpp

#include <pybind11/pybind11.h>

namespace py = pybind11;

// Forward declarations of binding initialization functions
void init_basic_bindings(py::module& m);
void init_parser_bindings(py::module& m);
void init_structure_bindings(py::module& m);
void init_forcefield_bindings(py::module& m);
void init_project_bindings(py::module& m);
void init_system_bindings(py::module& m);

PYBIND11_MODULE(pygcmc, m) {
    m.doc() = "Python bindings for GCMC library";

    // Initialize all bindings
    init_basic_bindings(m);
    init_parser_bindings(m);
    init_structure_bindings(m);
    init_forcefield_bindings(m);
    init_project_bindings(m);
    init_system_bindings(m);
}
