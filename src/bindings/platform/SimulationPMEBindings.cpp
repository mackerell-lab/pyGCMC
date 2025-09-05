// src/bindings/simulation/SimulationPMEBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Use new EnergyAPI instead of simulation.hpp
#include "../../platform/cpu/energy/EnergyAPI.hpp"
#include "../../model/ModelModule.hpp"
using namespace pygcmc;

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace platform {

void init_pme_bindings(py::module& m) {
    // PME bindings
    m.def("setPMEParameters",
        [](float alpha, const std::vector<int>& meshSize, int splineOrder, float tolerance) {
            if (meshSize.size() != 3) {
                throw std::runtime_error("meshSize must have exactly three elements");
            }
            int meshSize_array[3] = { meshSize[0], meshSize[1], meshSize[2] };
            ::pygcmc::platform::cpu::energy::setPMEParameters(alpha, meshSize_array, splineOrder, tolerance);
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
                ::pygcmc::platform::cpu::energy::initializePMEParameters(cutoff, box_array, alpha, nullptr, splineOrder, tolerance);
            } else {
                // If meshSize is provided, check size and convert to array
                if (meshSize.size() != 3) {
                    throw std::runtime_error("meshSize must have exactly three elements");
                }
                int meshSize_array[3] = {meshSize[0], meshSize[1], meshSize[2]};
                ::pygcmc::platform::cpu::energy::initializePMEParameters(cutoff, box_array, alpha, meshSize_array, splineOrder, tolerance);
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
            ::pygcmc::platform::cpu::energy::computeSystemEnergyPME(state);
            
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
            ::pygcmc::platform::cpu::energy::computeMovementEnergyPME(state);
            
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
        
    // Fixed versions that correct LJ double-counting
    m.def("computeSystemEnergyPMEFixed", 
        [](::pygcmc::model::MCState& state) {
            // Call C++ function to calculate energy
            ::pygcmc::platform::cpu::energy::computeSystemEnergyPME(state);
            
            // Convert from C++ struct to Python dictionary
            py::dict pme_dict;
            pme_dict["real_space"] = state.ewald_energy.real_space;
            pme_dict["reciprocal"] = state.ewald_energy.reciprocal;
            pme_dict["self"] = state.ewald_energy.self;
            
            // Calculate total electrostatic energy
            double electrostatic_total = state.ewald_energy.real_space + 
                                         state.ewald_energy.reciprocal + 
                                         state.ewald_energy.self;
            
            // Accumulate VDW energy from residues and fix double-counting
            double vdw = 0.0;
            for(const auto& res : state.residues) {
                if(res.active) {
                    vdw += res.energy_vdw;
                }
            }
            // Fix double-counting: divide by 2
            vdw /= 2.0;
            
            // Also update residue energies to be corrected
            for(auto& res : state.residues) {
                if(res.active) {
                    res.energy_vdw /= 2.0;
                }
            }
            
            // Correctly calculate and save total energy
            double total = electrostatic_total + vdw;
            pme_dict["total"] = total;
            
            // Return tuple: (electrostatic_total, vdw_energy, pme_dict)
            return py::make_tuple(electrostatic_total, vdw, pme_dict);
        },
        "Calculate system energy using PME with LJ double-counting fix");
        
    m.def("computeSystemEnergyCutoffFixed", 
        [](::pygcmc::model::MCState& state) {
            // Call the standard cutoff computation
            ::pygcmc::platform::cpu::energy::computeSystemEnergyCutoff(state);
            
            // Fix double-counting by dividing residue energies by 2
            for(auto& res : state.residues) {
                if(res.active) {
                    res.energy_vdw /= 2.0;
                    res.energy_elec /= 2.0;
                }
            }
            
            // Calculate corrected totals
            double total_vdw = 0.0;
            double total_elec = 0.0;
            for(const auto& res : state.residues) {
                if(res.active) {
                    total_vdw += res.energy_vdw;
                    total_elec += res.energy_elec;
                }
            }
            
            // Return tuple: (total_elec, total_vdw, total_energy)
            return py::make_tuple(total_elec, total_vdw, total_elec + total_vdw);
        },
        "Calculate cutoff energy with LJ double-counting fix");
        
    m.def("computeMovementEnergyPMEFixed", 
        [](::pygcmc::model::MCState& state) {
            // Call C++ function to calculate energy
            ::pygcmc::platform::cpu::energy::computeMovementEnergyPME(state);
            
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
            // Fix double-counting
            vdw /= 2.0;
            
            // Also update movement residue energies
            for(const auto& movementInfo : state.movementResidues) {
                for(int i = movementInfo.startIndex;
                    i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                    if(state.residues[i].active) {
                        state.residues[i].energy_vdw /= 2.0;
                    }
                }
            }
            
            // Correctly calculate and save total energy
            double total = electrostatic_total + vdw;
            pme_dict["total"] = total;
            
            // Return tuple: (electrostatic_total, vdw_energy, pme_dict)
            return py::make_tuple(electrostatic_total, vdw, pme_dict);
        },
        "Calculate movement residue energy using PME with LJ double-counting fix");
        
    // Complete energy calculation functions
    m.def("computeSystemEnergyPMEComplete", 
        [](::pygcmc::model::MCState& state) {
            ::pygcmc::platform::cpu::energy::computeSystemEnergyPMEComplete(state);
            
            // Extract energy components  
            double elec = state.ewald_energy.real_space + state.ewald_energy.reciprocal + state.ewald_energy.self;
            double vdw = 0.0;
            
            // Sum VdW energy from residues (no division by 2 needed for Complete)
            for(auto& res : state.residues) {
                if(res.active) {
                    vdw += res.energy_vdw;
                }
            }
            
            double total = elec + vdw;
            
            return std::make_tuple(elec, vdw, total);
        },
        py::arg("state"),
        "Compute complete system energy using PME with all interactions including intramolecular");
        
    m.def("computeSystemEnergyCutoffComplete",
        [](::pygcmc::model::MCState& state) {
            ::pygcmc::platform::cpu::energy::computeSystemEnergyCutoffComplete(state);
            
            // Extract energy components from C++ struct
            double elec = state.ewald_energy.real_space + state.ewald_energy.self;
            double vdw = 0.0;
            
            // Sum VdW energy from residues (no division by 2 needed for Complete)
            for(auto& res : state.residues) {
                if(res.active) {
                    vdw += res.energy_vdw;
                }
            }
            
            double total = elec + vdw;
            
            return std::make_tuple(elec, vdw, total);
        },
        py::arg("state"),
        "Compute complete system energy using cutoff with all interactions including intramolecular");
}

} // namespace platform
} // namespace bindings
} // namespace pygcmc