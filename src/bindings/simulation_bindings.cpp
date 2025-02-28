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
            // 调用C++函数计算能量
            simulation::Simulation::computeSystemEnergyEwald(state);
            
            // 从C++结构体转换为Python字典
            py::dict ewald_dict;
            ewald_dict["real_space"] = state.ewald_energy.real_space;
            ewald_dict["reciprocal"] = state.ewald_energy.reciprocal;
            ewald_dict["self"] = state.ewald_energy.self;
            
            // 计算静电总能量
            double electrostatic_total = state.ewald_energy.real_space + 
                                         state.ewald_energy.reciprocal + 
                                         state.ewald_energy.self;
            
            // 从残基中累计VDW能量
            double vdw = 0.0;
            for(const auto& res : state.residues) {
                if(res.active) {
                    vdw += res.energy_vdw;
                }
            }
            
            // 正确计算总能量并保存
            double total = electrostatic_total + vdw;
            ewald_dict["total"] = total;
            
            // 返回元组：(静电总能量, 范德华能量, ewald字典)
            return py::make_tuple(electrostatic_total, vdw, ewald_dict);
        },
        "Calculate system energy using Ewald summation");
          
    m.def("computeMovementEnergyEwald", 
        [](model::MCState& state) {
            // 调用C++函数计算能量
            simulation::Simulation::computeMovementEnergyEwald(state);
            
            // 从C++结构体转换为Python字典
            py::dict ewald_dict;
            ewald_dict["real_space"] = state.ewald_energy.real_space;
            ewald_dict["reciprocal"] = state.ewald_energy.reciprocal;
            ewald_dict["self"] = state.ewald_energy.self;
            
            // 计算静电总能量
            double electrostatic_total = state.ewald_energy.real_space + 
                                         state.ewald_energy.reciprocal + 
                                         state.ewald_energy.self;
            
            // 只累计movement残基的VDW能量
            double vdw = 0.0;
            for(const auto& movementInfo : state.movementResidues) {
                for(int i = movementInfo.startIndex;
                    i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                    if(state.residues[i].active) {
                        vdw += state.residues[i].energy_vdw;
                    }
                }
            }
            
            // 正确计算总能量并保存
            double total = electrostatic_total + vdw;
            ewald_dict["total"] = total;
            
            // 返回元组：(静电总能量, 范德华能量, ewald字典)
            return py::make_tuple(electrostatic_total, vdw, ewald_dict);
        },
        "Calculate movement residue energy using Ewald summation");
}

} // namespace bindings
} // namespace pygcmc
