#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../simulation/simulation.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {

void init_simulation_bindings(py::module& m) {
    py::enum_<platform::LogLevel>(m, "PlatformLogLevel")
        .value("DEBUG", platform::LogLevel::DEBUG)
        .value("INFO", platform::LogLevel::INFO)
        .value("WARNING", platform::LogLevel::WARNING)
        .value("ERROR", platform::LogLevel::ERROR);
    
    m.def("set_platform_verbose", &platform::set_verbose, 
         "Set verbose mode for platform logging");
    m.def("set_platform_log_level", &platform::set_log_level,
         "Set the minimum log level for platform");
    m.def("set_platform_debug_mode", &simulation::set_debug_mode,
         "Enable or disable debug mode for platform logging");

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
    
    m.def("computeSystemVdwEnergyCutoff", &simulation::Simulation::computeSystemVdwEnergyCutoff,
          "Calculate VDW energies for the full system with distance cutoff");
          
    m.def("setEnergyDebugOutput", &simulation::Simulation::setEnergyDebugOutput,
          "Enable or disable debug output for energy calculations");
    
    // CHARMM switching function related bindings have been removed
    // Please use the set_switching_function and calculate_switching_function methods in the MonteCarloSystem class

    // Ewald parameters are stored as static variables in the implementation
    m.def("setEwaldParameters",
        [](float alpha, const std::vector<int>& kmax, float tolerance) {
            if (kmax.size() != 3) {
                throw std::runtime_error("kmax must have exactly three elements");
            }
            int kmax_array[3] = { kmax[0], kmax[1], kmax[2] };
            // Implementation needs to be modified to make setEwaldParameters a static method
            simulation::Simulation::setEwaldParameters(alpha, kmax_array, tolerance);
        },
        "Set parameters for Ewald summation",
        py::arg("alpha"),
        py::arg("kmax"),
        py::arg("tolerance") = 1e-5f);
    
    m.def("initializeEwaldParameters",
        [](float cutoff, const std::vector<float>& box, float alpha, float tolerance) {
            if (box.size() != 3) {
                throw std::runtime_error("box must have exactly three elements");
            }
            // Convert to array and call function
            float box_array[3] = {box[0], box[1], box[2]};
            simulation::Simulation::initializeEwaldParameters(cutoff, box_array, alpha, tolerance);
        },
        "Initialize Ewald parameters with automatic optimization",
        py::arg("cutoff"),
        py::arg("box"),
        py::arg("alpha") = 0.0f,
        py::arg("tolerance") = 1e-5f);
          
    m.def("computeSystemEnergyEwald", 
        [](model::MCState& state) {
            // Call C++ function to calculate energy
            simulation::Simulation::computeSystemEnergyEwald(state);
            
            // Convert from C++ struct to Python dictionary
            py::dict ewald_dict;
            ewald_dict["real_space"] = state.ewald_energy.real_space;
            ewald_dict["reciprocal"] = state.ewald_energy.reciprocal;
            ewald_dict["self"] = state.ewald_energy.self;
            
            // Calculate total electrostatic energy
            double electrostatic_total = state.ewald_energy.real_space + 
                                         state.ewald_energy.reciprocal + 
                                         state.ewald_energy.self;
            
            // Accumulate VDW energy from residues
            double vdw = 0.0;
            for(const auto& res : state.residues) {
                if(res.active) {
                    vdw += res.energy_vdw;
                }
            }
            
            // Correctly calculate and save total energy
            double total = electrostatic_total + vdw;
            ewald_dict["total"] = total;
            
            // Return tuple: (electrostatic_total, vdw_energy, ewald_dict)
            return py::make_tuple(electrostatic_total, vdw, ewald_dict);
        },
        "Calculate system energy using Ewald summation");
          
    m.def("computeMovementEnergyEwald", 
        [](model::MCState& state) {
            // Call C++ function to calculate energy
            simulation::Simulation::computeMovementEnergyEwald(state);
            
            // Convert from C++ struct to Python dictionary
            py::dict ewald_dict;
            ewald_dict["real_space"] = state.ewald_energy.real_space;
            ewald_dict["reciprocal"] = state.ewald_energy.reciprocal;
            ewald_dict["self"] = state.ewald_energy.self;
            
            // Calculate total electrostatic energy
            double electrostatic_total = state.ewald_energy.real_space + 
                                         state.ewald_energy.reciprocal + 
                                         state.ewald_energy.self;
            
            // Only accumulate VDW energy from movement residues
            double vdw = 0.0;
            for(const auto& movementInfo : state.movementResidues) {
                for(int i = movementInfo.startIndex;
                    i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                    if(state.residues[i].active) {
                        vdw += state.residues[i].energy_vdw;
                    }
                }
            }
            
            // Correctly calculate and save total energy
            double total = electrostatic_total + vdw;
            ewald_dict["total"] = total;
            
            // Return tuple: (electrostatic_total, vdw_energy, ewald_dict)
            return py::make_tuple(electrostatic_total, vdw, ewald_dict);
        },
        "Calculate movement residue energy using Ewald summation");
    
    // PME bindings
    m.def("setPMEParameters",
        [](float alpha, const std::vector<int>& meshSize, int splineOrder, float tolerance) {
            if (meshSize.size() != 3) {
                throw std::runtime_error("meshSize must have exactly three elements");
            }
            int meshSize_array[3] = { meshSize[0], meshSize[1], meshSize[2] };
            simulation::Simulation::setPMEParameters(alpha, meshSize_array, splineOrder, tolerance);
        },
        "Set parameters for Particle Mesh Ewald summation",
        py::arg("alpha"),
        py::arg("meshSize"),
        py::arg("splineOrder") = 4,
        py::arg("tolerance") = 1e-5f);
    
    m.def("initializePMEParameters",
        [](float cutoff, const std::vector<float>& box, float alpha, 
           const std::vector<int>& meshSize, int splineOrder, float tolerance) {
            if (box.size() != 3) {
                throw std::runtime_error("box must have exactly three elements");
            }
            
            // Convert box to array
            float box_array[3] = {box[0], box[1], box[2]};
            
            // Handle optional meshSize
            if (meshSize.empty()) {
                // If meshSize is not provided, pass nullptr
                simulation::Simulation::initializePMEParameters(cutoff, box_array, alpha, nullptr, splineOrder, tolerance);
            } else {
                // If meshSize is provided, check size and convert to array
                if (meshSize.size() != 3) {
                    throw std::runtime_error("meshSize must have exactly three elements");
                }
                int meshSize_array[3] = {meshSize[0], meshSize[1], meshSize[2]};
                simulation::Simulation::initializePMEParameters(cutoff, box_array, alpha, meshSize_array, splineOrder, tolerance);
            }
        },
        "Initialize PME parameters with automatic optimization",
        py::arg("cutoff"),
        py::arg("box"),
        py::arg("alpha") = 0.0f,
        py::arg("meshSize") = std::vector<int>(),
        py::arg("splineOrder") = 4,
        py::arg("tolerance") = 1e-5f);
          
    m.def("computeSystemEnergyPME", 
        [](model::MCState& state) {
            // Call C++ function to calculate energy
            simulation::Simulation::computeSystemEnergyPME(state);
            
            // Convert from C++ struct to Python dictionary
            py::dict pme_dict;
            pme_dict["real_space"] = state.ewald_energy.real_space;
            pme_dict["reciprocal"] = state.ewald_energy.reciprocal;
            pme_dict["self"] = state.ewald_energy.self;
            
            // Calculate total electrostatic energy
            double electrostatic_total = state.ewald_energy.real_space + 
                                         state.ewald_energy.reciprocal + 
                                         state.ewald_energy.self;
            
            // Accumulate VDW energy from residues
            double vdw = 0.0;
            for(const auto& res : state.residues) {
                if(res.active) {
                    vdw += res.energy_vdw;
                }
            }
            
            // Correctly calculate and save total energy
            double total = electrostatic_total + vdw;
            pme_dict["total"] = total;
            
            // Return tuple: (electrostatic_total, vdw_energy, pme_dict)
            return py::make_tuple(electrostatic_total, vdw, pme_dict);
        },
        "Calculate system energy using Particle Mesh Ewald summation");
          
    m.def("computeMovementEnergyPME", 
        [](model::MCState& state) {
            // Call C++ function to calculate energy
            simulation::Simulation::computeMovementEnergyPME(state);
            
            // Convert from C++ struct to Python dictionary
            py::dict pme_dict;
            pme_dict["real_space"] = state.ewald_energy.real_space;
            pme_dict["reciprocal"] = state.ewald_energy.reciprocal;
            pme_dict["self"] = state.ewald_energy.self;
            
            // Calculate total electrostatic energy
            double electrostatic_total = state.ewald_energy.real_space + 
                                         state.ewald_energy.reciprocal + 
                                         state.ewald_energy.self;
            
            // Only accumulate VDW energy from movement residues
            double vdw = 0.0;
            for(const auto& movementInfo : state.movementResidues) {
                for(int i = movementInfo.startIndex;
                    i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                    if(state.residues[i].active) {
                        vdw += state.residues[i].energy_vdw;
                    }
                }
            }
            
            // Correctly calculate and save total energy
            double total = electrostatic_total + vdw;
            pme_dict["total"] = total;
            
            // Return tuple: (electrostatic_total, vdw_energy, pme_dict)
            return py::make_tuple(electrostatic_total, vdw, pme_dict);
        },
        "Calculate movement residue energy using Particle Mesh Ewald summation");
        
    // PGP bindings
    m.def("setPGPParameters",
        [](float alpha, const std::vector<int>& meshSize, float potential_cutoff, const std::vector<int>& potentialGridSize, 
           int splineOrder, float tolerance) {
            if (meshSize.size() != 3) {
                throw std::runtime_error("meshSize must have exactly three elements");
            }
            if (potentialGridSize.size() != 3) {
                throw std::runtime_error("potentialGridSize must have exactly three elements");
            }
            
            int meshSize_array[3] = { meshSize[0], meshSize[1], meshSize[2] };
            int potentialGridSize_array[3] = { potentialGridSize[0], potentialGridSize[1], potentialGridSize[2] };
            
            simulation::Simulation::setPGPParameters(alpha, meshSize_array, potential_cutoff, 
                                                      potentialGridSize_array, splineOrder, tolerance);
        },
        "Set parameters for Precomputed Grid-Potential PME summation",
        py::arg("alpha"),
        py::arg("meshSize"),
        py::arg("potential_cutoff"),
        py::arg("potentialGridSize"),
        py::arg("splineOrder") = 4,
        py::arg("tolerance") = 1e-5f,
        R"docstring(
        设置PGP-PME (Precomputed Grid-Potential Particle Mesh Ewald)算法参数

        PGP-PME是一种针对蒙特卡洛模拟优化的长程静电相互作用计算方法。它通过预计算
        系统中固定部分的静电势网格，大大加速了能量评估过程。

        参数:
            alpha (float): Ewald分离参数，控制实空间和倒空间计算的平衡
            meshSize (list[int]): PME网格尺寸 [nx, ny, nz]
            potential_cutoff (float): 电势计算的截断距离
            potentialGridSize (list[int]): 预计算电势网格尺寸 [nx, ny, nz]
            splineOrder (int, optional): B样条插值阶数，默认为4
            tolerance (float, optional): 精度容限，默认为1e-5
        )docstring");
        
    // 新增的PGP核心函数绑定
    m.def("precomputeGridPotential",
        [](model::MCState& state, bool fixed_only) {
            simulation::Simulation::precomputeGridPotential(state, fixed_only);
        },
        "Precompute the electrostatic grid potential for fixed parts of the system",
        py::arg("state"),
        py::arg("fixed_only") = true,
        R"docstring(
        预计算系统中固定部分的电势网格 (Precomputed Grid-Potential Particle Mesh Ewald)

        这是PGP-PME算法的核心函数之一，负责计算并存储系统中固定部分的静电势场。
        该预计算步骤只需在系统固定部分发生变化时执行一次，显著提高蒙特卡洛模拟效率。

        参数:
            state (MCState): 系统状态，包含原子坐标、电荷和盒子信息
            fixed_only (bool, optional): 是否只计算固定部分，默认为True
        )docstring");
        
    m.def("interpolateMoleculeEnergy",
        [](model::MCState& state) {
            double energy = 0.0;
            simulation::Simulation::interpolateMoleculeEnergy(state, energy);
            return energy;
        },
        "Calculate molecule energy by interpolating from the precomputed grid potential",
        py::arg("state"),
        R"docstring(
        通过插值计算移动分子的能量 (Precomputed Grid-Potential Particle Mesh Ewald)

        这是PGP-PME算法的另一个核心函数，通过B样条插值从预计算的电势网格中
        快速评估移动分子的能量，避免了直接计算分子间相互作用。

        参数:
            state (MCState): 系统状态，包含移动分子信息
        
        返回:
            float: 计算得到的能量值
        )docstring");
        
    // 添加新的函数绑定：calculateMoleculeEnergy
    m.def("calculateMoleculeEnergy",
        [](model::MCState& state) {
            return simulation::Simulation::calculateMoleculeEnergy(state);
        },
        "Calculate molecule energy by interpolating from the precomputed grid potential (alternative function)",
        py::arg("state"),
        R"docstring(
        通过插值计算移动分子的能量 (Precomputed Grid-Potential Particle Mesh Ewald)

        这是与interpolateMoleculeEnergy等价的函数，提供更直观的命名。
        通过B样条插值从预计算的电势网格中获取移动分子的能量值。
        
        参数:
            state (MCState): 系统状态，包含移动分子信息
        
        返回:
            float: 计算得到的能量值
        )docstring");
}

} // namespace bindings
} // namespace pygcmc
