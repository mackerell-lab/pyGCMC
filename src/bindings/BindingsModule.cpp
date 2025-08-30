// src/bindings/BindingsModule.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "platform/cpu/energy/common/MemoryGlobalCleanup.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {

// Forward declarations for Simulation binding init functions from separate files
namespace simulation {
void init_basic_bindings(py::module& m);
void init_ewald_bindings(py::module& m);
void init_pme_bindings(py::module& m);
void init_pgp_bindings(py::module& m);
void init_drude_bindings(py::module& m);
void init_movement_bindings(py::module& m);
void init_gcmc_bindings(py::module& m);

void init_simulation_bindings(py::module& m) {
    // Initialize all Simulation binding groups
    init_basic_bindings(m);
    init_ewald_bindings(m);
    init_pme_bindings(m);
    init_pgp_bindings(m);
    init_drude_bindings(m);
    init_movement_bindings(m);
    init_gcmc_bindings(m);
}
}

// Forward declarations for System binding init functions from separate files
namespace system {
void init_common_bindings(py::module& m, py::module& system_module);
void init_molecular_bindings(py::module& m, py::module& system_module);
void init_montecarlo_bindings(py::module& m, py::module& system_module);

void init_system(py::module& m) {
    // Create system submodule
    auto system_module = m.def_submodule("system", "System management classes");
    
    // Initialize all System binding groups
    init_common_bindings(m, system_module);
    init_molecular_bindings(m, system_module);
    init_montecarlo_bindings(m, system_module);
}
}

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

// Forward declarations for Platform binding init functions
namespace platform {
void init_energy_bindings(py::module& m);
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

// Cleanup function to be called before Python exits
static void cleanup_global_state() {
    pygcmc::platform::cpu::cleanupAllGlobalState();
}

// Main pybind11 module definition
PYBIND11_MODULE(pygcmc, m) {
    m.doc() = "Python bindings for GCMC simulation library";
    
    // Initialize bindings - io, model, and system are refactored
    pygcmc::bindings::io::init_io_bindings(m);
    pygcmc::bindings::model::init_model(m);
    pygcmc::bindings::system::init_system(m);
    pygcmc::bindings::simulation::init_simulation_bindings(m);
    pygcmc::bindings::platform::init_energy_bindings(m);
    
    // Register cleanup function - can be called manually if needed
    // NOTE: We do NOT automatically register with atexit to avoid 
    // pybind11 deallocation issues during Python shutdown
    m.def("_cleanup", &cleanup_global_state, 
          "Internal cleanup function - call manually before exit if needed",
          py::call_guard<py::gil_scoped_release>());
    
    // Export OpenMP capability flag
#ifdef PYGCMC_USE_OPENMP
    m.attr("PYGCMC_USE_OPENMP") = py::bool_(true);
#else
    m.attr("PYGCMC_USE_OPENMP") = py::bool_(false);
#endif
}