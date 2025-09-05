// src/bindings/simulation/SimulationPGPBindings.cpp

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

void init_pgp_bindings(py::module& m) {
    // PGP bindings
    m.def("resetPGPState",
        []() {
            ::pygcmc::platform::cpu::energy::resetPGPState();
        },
        "Reset PGP global state to fix memory corruption issues. "
        "Call this between tests or when reinitializing PGP parameters.");
        
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
            
            ::pygcmc::platform::cpu::energy::setPGPParameters(alpha, meshSize_array, potential_cutoff, 
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
        [](::pygcmc::model::MCState& state, bool fixed_only) {
            ::pygcmc::platform::cpu::energy::precomputeGridPotential(state, fixed_only);
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
        [](::pygcmc::model::MCState& state) {
            double energy = 0.0;
            ::pygcmc::platform::cpu::energy::interpolateMoleculeEnergy(state, energy);
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
        [](::pygcmc::model::MCState& state) {
            return ::pygcmc::platform::cpu::energy::calculateMoleculeEnergy(state);
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
        [](::pygcmc::model::MCState& state) {
            // Call fixed C++ function to calculate energy
            ::pygcmc::platform::cpu::energy::computeSystemEnergyPGPFixed(state);
            
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
        [](::pygcmc::model::MCState& state) {
            // Call fixed C++ function to calculate energy
            ::pygcmc::platform::cpu::energy::computeMovementEnergyPGPFixed(state);
            
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
        
    // PGP Complete bindings
    m.def("computeSystemEnergyPGPComplete", 
        [](::pygcmc::model::MCState& state) {
            // Call C++ function
            ::pygcmc::platform::cpu::energy::computeSystemEnergyPGPComplete(state);
            
            // Calculate total electrostatic energy
            double electrostatic_total = state.ewald_energy.real_space + 
                                       state.ewald_energy.reciprocal + 
                                       state.ewald_energy.self;
            
            // Accumulate VDW energy from all residues
            double vdw = 0.0;
            for (const auto& res : state.residues) {
                if (res.active) {
                    vdw += res.energy_vdw;
                }
            }
            
            // Return tuple: (electrostatic_total, vdw_energy, total_energy)
            return py::make_tuple(electrostatic_total, vdw, electrostatic_total + vdw);
        },
        "Calculate complete system energy using PGP with intramolecular LJ interactions");
        
    m.def("computeMovementEnergyPGPComplete", 
        [](::pygcmc::model::MCState& state) {
            // Call C++ function
            ::pygcmc::platform::cpu::energy::computeMovementEnergyPGPComplete(state);
            
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
            
            // Save total energy
            double total = electrostatic_total + vdw;
            pgp_dict["total"] = total;
            
            // Return tuple: (electrostatic_total, vdw_energy, pgp_dict)
            return py::make_tuple(electrostatic_total, vdw, pgp_dict);
        },
        "Calculate movement energy using PGP Complete with LJ interactions");
        
    m.def("computeMovementEnergyPGPCompleteCorrect", 
        [](::pygcmc::model::MCState& state) {
            // Call existing function (the "Correct" version doesn't exist yet)
            ::pygcmc::platform::cpu::energy::computeMovementEnergyPGPComplete(state);
            
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
            if(state.movementResidues.empty()) {
                // If no movement residues specified, sum VDW from all non-fixed residues
                for(int i = 0; i < state.activeResidueCount; i++) {
                    if(state.residues[i].active && !state.residues[i].fixed) {
                        vdw += state.residues[i].energy_vdw;
                    }
                }
            } else {
                for(const auto& movementInfo : state.movementResidues) {
                    for(int i = movementInfo.startIndex;
                        i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                        if(state.residues[i].active) {
                            vdw += state.residues[i].energy_vdw;
                        }
                    }
                }
            }
            
            // Save total energy
            double total = electrostatic_total + vdw;
            pgp_dict["total"] = total;
            
            // Return tuple: (electrostatic_total, vdw_energy, pgp_dict)
            return py::make_tuple(electrostatic_total, vdw, pgp_dict);
        },
        "Corrected calculation of movement energy using PGP Complete");
}

} // namespace platform
} // namespace bindings
} // namespace pygcmc