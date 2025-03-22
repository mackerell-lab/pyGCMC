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
        [](float alpha, const std::vector<int>& meshSize, float pair_cutoff, const std::vector<int>& pairGridSize, 
           int splineOrder, float tolerance) {
            if (meshSize.size() != 3) {
                throw std::runtime_error("meshSize must have exactly three elements");
            }
            if (pairGridSize.size() != 3) {
                throw std::runtime_error("pairGridSize must have exactly three elements");
            }
            
            int meshSize_array[3] = { meshSize[0], meshSize[1], meshSize[2] };
            int pairGridSize_array[3] = { pairGridSize[0], pairGridSize[1], pairGridSize[2] };
            
            simulation::Simulation::setPGPParameters(alpha, meshSize_array, pair_cutoff, 
                                                      pairGridSize_array, splineOrder, tolerance);
        },
        "Set parameters for Pair-Grid PME summation",
        py::arg("alpha"),
        py::arg("meshSize"),
        py::arg("pair_cutoff"),
        py::arg("pairGridSize"),
        py::arg("splineOrder") = 4,
        py::arg("tolerance") = 1e-5f);
    
    m.def("initializePGPParameters",
        [](float cutoff, float pair_cutoff, const std::vector<float>& box, float alpha, 
           const std::vector<int>& meshSize, const std::vector<int>& pairGridSize, 
           int splineOrder, float tolerance) {
            if (box.size() != 3) {
                throw std::runtime_error("box must have exactly three elements");
            }
            
            // Convert box to array
            float box_array[3] = {box[0], box[1], box[2]};
            
            // Handle optional meshSize and pairGridSize
            const int* meshSize_ptr = nullptr;
            const int* pairGridSize_ptr = nullptr;
            
            // If meshSize is provided, check size and convert to array
            int meshSize_array[3] = {0, 0, 0};
            if (!meshSize.empty()) {
                if (meshSize.size() != 3) {
                    throw std::runtime_error("meshSize must have exactly three elements");
                }
                meshSize_array[0] = meshSize[0];
                meshSize_array[1] = meshSize[1];
                meshSize_array[2] = meshSize[2];
                meshSize_ptr = meshSize_array;
            }
            
            // If pairGridSize is provided, check size and convert to array
            int pairGridSize_array[3] = {0, 0, 0};
            if (!pairGridSize.empty()) {
                if (pairGridSize.size() != 3) {
                    throw std::runtime_error("pairGridSize must have exactly three elements");
                }
                pairGridSize_array[0] = pairGridSize[0];
                pairGridSize_array[1] = pairGridSize[1];
                pairGridSize_array[2] = pairGridSize[2];
                pairGridSize_ptr = pairGridSize_array;
            }
            
            simulation::Simulation::initializePGPParameters(cutoff, pair_cutoff, box_array, alpha, 
                                                             meshSize_ptr, pairGridSize_ptr, 
                                                             splineOrder, tolerance);
        },
        "Initialize PGP parameters with automatic optimization",
        py::arg("cutoff"),
        py::arg("pair_cutoff"),
        py::arg("box"),
        py::arg("alpha") = 0.0f,
        py::arg("meshSize") = std::vector<int>(),
        py::arg("pairGridSize") = std::vector<int>(),
        py::arg("splineOrder") = 4,
        py::arg("tolerance") = 1e-5f);
          
    m.def("computeSystemEnergyPGP", 
        [](model::MCState& state) {
            // Call C++ function to calculate energy
            simulation::Simulation::computeSystemEnergyPGP(state);
            
            // Convert from C++ struct to Python dictionary
            py::dict pgp_dict;
            pgp_dict["real_space"] = state.ewald_energy.real_space;
            pgp_dict["reciprocal"] = state.ewald_energy.reciprocal;
            pgp_dict["self"] = state.ewald_energy.self;
            
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
            pgp_dict["total"] = total;
            
            // Return tuple: (electrostatic_total, vdw_energy, pgp_dict)
            return py::make_tuple(electrostatic_total, vdw, pgp_dict);
        },
        "Calculate system energy using Pair-Grid PME summation");
          
    m.def("computeMovementEnergyPGP", 
        [](model::MCState& state) {
            // Call C++ function to calculate energy
            simulation::Simulation::computeMovementEnergyPGP(state);
            
            // Convert from C++ struct to Python dictionary
            py::dict pgp_dict;
            pgp_dict["real_space"] = state.ewald_energy.real_space;
            pgp_dict["reciprocal"] = state.ewald_energy.reciprocal;
            pgp_dict["self"] = state.ewald_energy.self;
            
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
            pgp_dict["total"] = total;
            
            // Return tuple: (electrostatic_total, vdw_energy, pgp_dict)
            return py::make_tuple(electrostatic_total, vdw, pgp_dict);
        },
        "Calculate movement residue energy using Pair-Grid PME summation");
}

} // namespace bindings
} // namespace pygcmc
