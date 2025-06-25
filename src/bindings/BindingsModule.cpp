// src/bindings/BindingsModule.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace pygcmc {
namespace bindings {

// Forward declarations for binding init functions from other modules
void init_system(pybind11::module& m);  
void init_simulation_bindings(pybind11::module& m);

// Forward declarations for Model binding init functions from separate files
namespace model {
void init_structure_bindings(py::module& m, py::module& model_module);
void init_atom_bindings(py::module& m, py::module& model_module);
void init_residue_bindings(py::module& m, py::module& model_module);
void init_topology_bindings(py::module& m, py::module& model_module);
void init_forcefield_bindings(py::module& m, py::module& model_module);
void init_param_bindings(py::module& m, py::module& model_module);
void init_molecule_bindings(py::module& m, py::module& model_module);
void init_montecarlo_bindings(py::module& m, py::module& model_module);

void init_model(py::module& m) {
    // Create model submodule
    auto model_module = m.def_submodule("model", "Data model classes");
    
    // Initialize all Model binding groups
    init_structure_bindings(m, model_module);
    init_atom_bindings(m, model_module);
    init_residue_bindings(m, model_module);
    init_topology_bindings(m, model_module);
    init_forcefield_bindings(m, model_module);
    init_param_bindings(m, model_module);
    init_molecule_bindings(m, model_module);
    init_montecarlo_bindings(m, model_module);
}
}

// Forward declarations for IO binding init functions from separate files
namespace io {
void init_structure_bindings(py::module& m, py::module& io_module);
void init_topology_bindings(py::module& m, py::module& io_module);
void init_forcefield_bindings(py::module& m, py::module& io_module);
void init_parameters_bindings(py::module& m, py::module& io_module);

void init_io_bindings(py::module& m) {
    // Create io submodule
    auto io_module = m.def_submodule("io", "Input/Output operations");
    
    // Initialize all IO binding groups
    init_structure_bindings(m, io_module);
    init_topology_bindings(m, io_module);
    init_forcefield_bindings(m, io_module);
    init_parameters_bindings(m, io_module);
}

} // namespace io
} // namespace bindings
} // namespace pygcmc

// Main pybind11 module definition
PYBIND11_MODULE(pygcmc, m) {
    m.doc() = "Python bindings for GCMC simulation library";
    
    // Initialize bindings - io and model are refactored
    pygcmc::bindings::io::init_io_bindings(m);
    pygcmc::bindings::model::init_model(m);
    pygcmc::bindings::init_system(m);
    pygcmc::bindings::init_simulation_bindings(m);
}