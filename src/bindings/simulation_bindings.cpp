#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../simulation/simulation.hpp"
namespace py = pybind11;

namespace pygcmc {
namespace bindings {

void init_simulation_bindings(py::module& m) {
    // Bind the naive nonbonded energy calculation function
    m.def("computeNaiveNonbondedEnergy", &simulation::Simulation::computeNaiveNonbondedEnergy,
          "Compute naive nonbonded energy (fixed r=1.0, no PBC) between active movement residues and all other active residues",
          py::arg("state"));
}

} // namespace bindings
} // namespace pygcmc
