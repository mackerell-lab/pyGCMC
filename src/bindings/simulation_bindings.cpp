#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../simulation/simulation.hpp"
namespace py = pybind11;

namespace pygcmc {
namespace bindings {

void init_simulation_bindings(py::module& m) {
    // Bind the naive nonbonded energy calculation function
    m.def("computeNaiveNonbondedEnergy", &simulation::Simulation::computeNaiveNonbondedEnergy,
          "Compute and update nonbonded energies (vdw and elec) for all active movement residues",
          py::arg("state"));
          
    // 修改绑定，使用Simulation类的方法
    m.def("computeAllNonbondedEnergy", &simulation::Simulation::computeAllNonbondedEnergy,
          "Compute and update nonbonded energies (vdw and elec) for all active residues",
          py::arg("state"));
}

} // namespace bindings
} // namespace pygcmc
