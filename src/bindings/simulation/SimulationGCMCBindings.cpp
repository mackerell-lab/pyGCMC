// GCMCBindings.cpp - Simplified Python bindings for GCMC simulation

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../simulation/simulation.hpp"
#include "../../platform/cpu/energy/common/EnergyInterface.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace simulation {

void init_gcmc_bindings(py::module& m) {
    using namespace ::pygcmc::simulation;
    
    // GCMCSimulation class
    py::class_<GCMCSimulation>(m, "GCMCSimulation", "GCMC simulation class")
        .def(py::init<const GCMCSimulation::Config&>(), py::arg("config"))
        .def("initialize", &GCMCSimulation::initialize, "Initialize with MC state",
             py::arg("state"))
        .def("add_water", &GCMCSimulation::addWater, "Add water template")
        .def("run", &GCMCSimulation::run, "Run simulation")
        .def("run_steps", &GCMCSimulation::runSteps, "Run specific number of steps",
             py::arg("n_steps"))
        .def("get_results", &GCMCSimulation::getResults, "Get simulation results");
    
    // GCMCSimulation::Config
    py::class_<GCMCSimulation::Config>(m, "GCMCSimulationConfig", "Configuration for GCMCSimulation")
        .def(py::init<>())
        .def_readwrite("temperature", &GCMCSimulation::Config::temperature)
        .def_readwrite("equilibration_steps", &GCMCSimulation::Config::equilibrationSteps)
        .def_readwrite("production_steps", &GCMCSimulation::Config::productionSteps)
        .def_readwrite("chemical_potential", &GCMCSimulation::Config::chemicalPotential)
        .def_readwrite("use_cavity_bias", &GCMCSimulation::Config::useCavityBias)
        .def_readwrite("verbose", &GCMCSimulation::Config::verbose);
    
    // GCMCSimulation::Results
    py::class_<GCMCSimulation::Results>(m, "GCMCResults", "Results from GCMC simulation")
        .def_readonly("average_molecules", &GCMCSimulation::Results::averageMolecules)
        .def_readonly("average_energy", &GCMCSimulation::Results::averageEnergy)
        .def_readonly("acceptance_rate", &GCMCSimulation::Results::acceptanceRate);
    
    // Energy method enum
    py::enum_<platform::cpu::EnergyMethod>(m, "EnergyMethod")
        .value("DIRECT", platform::cpu::EnergyMethod::DIRECT)
        .value("EWALD", platform::cpu::EnergyMethod::EWALD)
        .value("PME", platform::cpu::EnergyMethod::PME);
    
    // Quick GCMC functions
    m.def("run_water_simulation", &GCMC::runWaterSimulation,
          "Run water GCMC simulation with default parameters",
          py::arg("state"),
          py::arg("temperature") = 300.0,
          py::arg("chemical_potential") = -15.7,
          py::arg("steps") = 100000);
}

} // namespace simulation
} // namespace bindings
} // namespace pygcmc