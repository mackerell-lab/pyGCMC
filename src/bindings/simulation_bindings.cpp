#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../simulation/simulation.hpp"
namespace py = pybind11;

namespace pygcmc {
namespace bindings {

void init_simulation_bindings(py::module& m) {
    m.def("computeMovementResiduesEnergy", &simulation::Simulation::computeMovementResiduesEnergy,
          "Calculate nonbonded energies for movement residues only");
          
    m.def("computeFullSystemEnergy", &simulation::Simulation::computeFullSystemEnergy,
          "Calculate nonbonded energies for the full system");
          
    m.def("computeFullSystemCutoffEnergy", &simulation::Simulation::computeFullSystemCutoffEnergy,
          "Calculate nonbonded energies for the full system with distance cutoff");
          
    m.def("computeFullSystemCutoffPBCEnergy", &simulation::Simulation::computeFullSystemCutoffPBCEnergy,
          "Calculate nonbonded energies for the full system with distance cutoff and periodic boundary conditions");
}

} // namespace bindings
} // namespace pygcmc
