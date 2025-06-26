// src/bindings/simulation/SimulationPMEBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../simulation/simulation.hpp"
#include "../../model/ModelModule.hpp"
using namespace pygcmc;

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace simulation {

void init_pme_bindings(py::module& m) {
    // PME bindings
    m.def("setPMEParameters",
        [](float alpha, const std::vector<int>& meshSize, int splineOrder, float tolerance) {
            if (meshSize.size() != 3) {
                throw std::runtime_error("meshSize must have exactly three elements");
            }
            int meshSize_array[3] = { meshSize[0], meshSize[1], meshSize[2] };
            ::pygcmc::simulation::Simulation::setPMEParameters(alpha, meshSize_array, splineOrder, tolerance);
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
                ::pygcmc::simulation::Simulation::initializePMEParameters(cutoff, box_array, alpha, nullptr, splineOrder, tolerance);
            } else {
                // If meshSize is provided, check size and convert to array
                if (meshSize.size() != 3) {
                    throw std::runtime_error("meshSize must have exactly three elements");
                }
                int meshSize_array[3] = {meshSize[0], meshSize[1], meshSize[2]};
                ::pygcmc::simulation::Simulation::initializePMEParameters(cutoff, box_array, alpha, meshSize_array, splineOrder, tolerance);
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
        [](::pygcmc::model::MCState& state) {
            // Call C++ function to calculate energy
            ::pygcmc::simulation::Simulation::computeSystemEnergyPME(state);
            
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
        [](::pygcmc::model::MCState& state) {
            // Call C++ function to calculate energy
            ::pygcmc::simulation::Simulation::computeMovementEnergyPME(state);
            
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
}

} // namespace simulation
} // namespace bindings
} // namespace pygcmc