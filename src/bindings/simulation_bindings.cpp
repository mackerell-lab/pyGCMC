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
        Set PGP-PME (Precomputed Grid-Potential Particle Mesh Ewald) algorithm parameters

        PGP-PME is a long-range electrostatic interaction calculation method optimized for Monte Carlo simulations. It
        greatly accelerates the energy evaluation process by precomputing the electrostatic potential grid of the fixed parts of the system.

        Parameters:
            alpha (float): Ewald separation parameter, controls the balance between real-space and reciprocal-space calculations
            meshSize (list[int]): PME grid size [nx, ny, nz]
            potential_cutoff (float): Cutoff distance for potential calculation
            potentialGridSize (list[int]): Precomputed potential grid size [nx, ny, nz]
            splineOrder (int, optional): B-spline interpolation order, default is 4
            tolerance (float, optional): Precision tolerance, default is 1e-5
        )docstring");
        
    // Added core PGP function bindings
    m.def("precomputeGridPotential",
        [](model::MCState& state, bool fixed_only) {
            simulation::Simulation::precomputeGridPotential(state, fixed_only);
        },
        "Precompute the electrostatic potential grid for fixed parts of the system (Precomputed Grid-Potential Particle Mesh Ewald)",
        py::arg("state"),
        py::arg("fixed_only") = true,
        R"docstring(
        Precompute the electrostatic potential grid for fixed parts of the system (Precomputed Grid-Potential Particle Mesh Ewald)

        This is one of the core functions of the PGP-PME algorithm, responsible for calculating and storing the electrostatic potential field of the fixed parts of the system.
        This precomputation step only needs to be executed once when the fixed parts of the system change, significantly improving Monte Carlo simulation efficiency.

        Parameters:
            state (MCState): System state, containing atom coordinates, charges, and box information
            fixed_only (bool, optional): Whether to calculate only the fixed parts, default is True
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
        Calculate molecule energy by interpolating from the precomputed grid potential (Precomputed Grid-Potential Particle Mesh Ewald)

        This is another core function of the PGP-PME algorithm, using B-spline interpolation from the precomputed potential grid to
        quickly evaluate the energy of moving molecules, avoiding direct calculation of intermolecular interactions.

        Parameters:
            state (MCState): System state, containing information about moving molecules
        
        Returns:
            float: The calculated energy value
        )docstring");
        
    // Add new function binding: calculateMoleculeEnergy
    m.def("calculateMoleculeEnergy",
        [](model::MCState& state) {
            return simulation::Simulation::calculateMoleculeEnergy(state);
        },
        "Calculate molecule energy by interpolating from the precomputed grid potential (alternative function)",
        py::arg("state"),
        R"docstring(
        Calculate molecule energy by interpolating from the precomputed grid potential (Precomputed Grid-Potential Particle Mesh Ewald)

        This is an equivalent function to interpolateMoleculeEnergy, providing a more intuitive naming.
        It retrieves the energy value of moving molecules through B-spline interpolation from the precomputed potential grid.
        
        Parameters:
            state (MCState): System state, containing information about moving molecules
        
        Returns:
            float: The calculated energy value
        )docstring");

    m.def("computeSystemEnergyPGP", 
        [](model::MCState& state) {
            // Call C++ function to calculate energy
            simulation::Simulation::computeSystemEnergyPGP(state);
            
            // Convert from C++ struct to Python dictionary
            py::dict pgp_dict;
            pgp_dict["grid_energy"] = state.ewald_energy.reciprocal;  // grid energy is stored in reciprocal field
            pgp_dict["real_space"] = state.ewald_energy.real_space;
            pgp_dict["self"] = state.ewald_energy.self;
            
            // Calculate total electrostatic energy
            double electrostatic_total = state.ewald_energy.reciprocal + 
                                        state.ewald_energy.real_space + 
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
        "Calculate system energy using PGP-PME method");
          
    m.def("computeMovementEnergyPGP", 
        [](model::MCState& state) {
            // Call C++ function to calculate energy
            simulation::Simulation::computeMovementEnergyPGP(state);
            
            // Convert from C++ struct to Python dictionary
            py::dict pgp_dict;
            pgp_dict["grid_energy"] = state.ewald_energy.reciprocal;  // grid energy is stored in reciprocal field
            pgp_dict["real_space"] = state.ewald_energy.real_space;
            pgp_dict["self"] = state.ewald_energy.self;
            
            // Calculate total electrostatic energy
            double electrostatic_total = state.ewald_energy.reciprocal + 
                                        state.ewald_energy.real_space + 
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
        "Calculate movement residue energy using PGP-PME method");
}

} // namespace bindings
} // namespace pygcmc
