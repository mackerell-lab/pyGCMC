#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../simulation/simulation.hpp"
namespace py = pybind11;

namespace pygcmc {
namespace bindings {

void init_simulation_bindings(py::module& m) {
    m.def("computeMovementEnergy", &simulation::Simulation::computeMovementEnergy,
          "Calculate nonbonded energies for movement residues only");
          
    m.def("computeMovementEnergyCutoff", &simulation::Simulation::computeMovementEnergyCutoff,
          "Calculate nonbonded energies for movement residues only with distance cutoff");
          
    m.def("computeSystemEnergy", &simulation::Simulation::computeSystemEnergy,
          "Calculate nonbonded energies for the full system");
          
    m.def("computeSystemEnergyCutoff", &simulation::Simulation::computeSystemEnergyCutoff,
          "Calculate nonbonded energies for the full system with distance cutoff");
          
    m.def("computeSystemEnergyPBC", &simulation::Simulation::computeSystemEnergyPBC,
          "Calculate nonbonded energies for the full system with distance cutoff and periodic boundary conditions");
          
    m.def("setEnergyDebugOutput", &simulation::Simulation::setEnergyDebugOutput,
          "Enable or disable debug output for energy calculations");
}

} // namespace bindings
} // namespace pygcmc
