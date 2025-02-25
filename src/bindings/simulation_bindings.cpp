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
          "Calculate nonbonded energies for the full system with periodic boundary conditions");
          
    m.def("computeSystemEnergyPBCCutoff", &simulation::Simulation::computeSystemEnergyPBCCutoff,
          "Calculate nonbonded energies for the full system with periodic boundary conditions and cutoff");
          
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
          
    m.def("computeSystemEnergyEwald", 
        [](model::MCState& state) {
            simulation::Simulation::computeSystemEnergyEwald(state);
            
            // 计算总能量
            double electrostatic = 0.0;
            double vdw = 0.0;
            
            // 从ewald_energy结构体中获取能量
            electrostatic = state.ewald_energy.total;
            
            // 从残基中累计VDW能量
            for(const auto& res : state.residues) {
                if(res.active) {
                    vdw += res.energy_vdw;
                }
            }
            
            // 返回[静电能量, 范德华能量]的元组
            return py::make_tuple(electrostatic, vdw);
        },
        "Calculate system energy using Ewald summation");
          
    m.def("computeMovementEnergyEwald", &simulation::Simulation::computeMovementEnergyEwald,
          "Calculate movement residues energy using Ewald summation");
}

} // namespace bindings
} // namespace pygcmc
