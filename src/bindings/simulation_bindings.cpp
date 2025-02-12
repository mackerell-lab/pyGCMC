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

    // Ewald parameters are stored as static variables in the implementation
    m.def("setEwaldParameters",
        [](float alpha, const std::vector<int>& kmax, float tolerance) {
            if (kmax.size() != 3) {
                throw std::runtime_error("kmax must have exactly three elements");
            }
            int kmax_array[3] = { kmax[0], kmax[1], kmax[2] };
            // 这里需要修改实现，让setEwaldParameters成为静态方法
            simulation::Simulation::setEwaldParameters(alpha, kmax_array, tolerance);
        },
        "Set parameters for Ewald summation",
        py::arg("alpha"),
        py::arg("kmax"),
        py::arg("tolerance") = 1e-5f);
          
    m.def("computeSystemEnergyEwald", &simulation::Simulation::computeSystemEnergyEwald,
          "Calculate system energy using Ewald summation");
          
    m.def("computeMovementEnergyEwald", &simulation::Simulation::computeMovementEnergyEwald,
          "Calculate movement residues energy using Ewald summation");
}

} // namespace bindings
} // namespace pygcmc
