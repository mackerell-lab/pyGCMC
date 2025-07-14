// src/bindings/simulation/SimulationPGPBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <memory>
#include <mutex>

#include "../../simulation/simulation.hpp"
#include "../../model/ModelModule.hpp"
#include "../../platform/cpu/energy/pgp/PGPContext.hpp"
#include "../../platform/cpu/energy/pme/PMECore.hpp"

using namespace pygcmc;

namespace py = pybind11;

// Thread-safe global PGPContext instance for backward compatibility
static std::unique_ptr<::platform::cpu::energy::pgp::PGPContext> global_pgp_context;
static std::mutex pgp_mutex;
static bool pgp_initialized = false;

// Store PGP parameters for later initialization
struct PGPParams {
    float alpha = 2.5f;
    std::vector<int> meshSize = {32, 32, 32};
    float potential_cutoff = 10.0f;
    std::vector<int> potentialGridSize = {32, 32, 32};
    int splineOrder = 4;
    float tolerance = 1e-6f;
    bool set = false;
};
static PGPParams stored_pgp_params;

// Store PME parameters
struct PMEParams {
    float cutoff = 10.0f;
    float box[3] = {30.0f, 30.0f, 30.0f};
    float alpha = 2.5f;
    bool set = false;
};
static PMEParams stored_pme_params;

namespace pygcmc {
namespace bindings {
namespace simulation {

void init_pgp_bindings(py::module& m) {
    // PGP bindings
    m.def("resetPGPState",
        []() {
            std::lock_guard<std::mutex> lock(pgp_mutex);
            // Reset our global PGPContext and parameters
            global_pgp_context.reset();
            pgp_initialized = false;
            stored_pgp_params = PGPParams();
            stored_pme_params = PMEParams();
            
            // Also call the original reset function
            ::pygcmc::simulation::Simulation::resetPGPState();
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
            
            // Store parameters for later use with PGPContext
            std::lock_guard<std::mutex> lock(pgp_mutex);
            stored_pgp_params.alpha = alpha;
            stored_pgp_params.meshSize = meshSize;
            stored_pgp_params.potential_cutoff = potential_cutoff;
            stored_pgp_params.potentialGridSize = potentialGridSize;
            stored_pgp_params.splineOrder = splineOrder;
            stored_pgp_params.tolerance = tolerance;
            stored_pgp_params.set = true;
            
            // Also call the original function to set global state for PGP Complete
            int meshSize_array[3] = { meshSize[0], meshSize[1], meshSize[2] };
            int potentialGridSize_array[3] = { potentialGridSize[0], potentialGridSize[1], potentialGridSize[2] };
            ::pygcmc::simulation::Simulation::setPGPParameters(alpha, meshSize_array, potential_cutoff, 
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
        
    // Remove PGP-specific initializePMEParameters to avoid conflicts with PME binding
    // The PME binding version will set global PME parameters that we can read
    
    // Added core PGP function bindings
    m.def("precomputeGridPotential",
        [](::pygcmc::model::MCState& state, bool fixed_only [[maybe_unused]]) {
            std::lock_guard<std::mutex> lock(pgp_mutex);
            
            // Initialize PGPContext if not already done
            if (!pgp_initialized || !global_pgp_context) {
                // Debug output
                py::print("Debug: pgp_initialized =", pgp_initialized);
                py::print("Debug: global_pgp_context =", (global_pgp_context ? "set" : "null"));
                py::print("Debug: stored_pme_params.set =", stored_pme_params.set);
                py::print("Debug: stored_pme_params.cutoff =", stored_pme_params.cutoff);
                py::print("Debug: stored_pme_params.box =", 
                         "[", stored_pme_params.box[0], ",", 
                         stored_pme_params.box[1], ",", 
                         stored_pme_params.box[2], "]");
                py::print("Debug: stored_pgp_params.set =", stored_pgp_params.set);
                py::print("Debug: stored_pgp_params.alpha =", stored_pgp_params.alpha);
                
                // Always read from global PME parameters
                bool pme_initialized_globally = false;
                try {
                    // Check if global PME parameters are available
                    auto& pme = ::pygcmc::platform::cpu::pme_params;
                    if (pme.initialized && pme.meshSize[0] > 0) {
                        pme_initialized_globally = true;
                        // Always use global PME parameters
                        stored_pme_params.cutoff = pme.cutoff;
                        stored_pme_params.box[0] = static_cast<float>(pme.box[0]);
                        stored_pme_params.box[1] = static_cast<float>(pme.box[1]);
                        stored_pme_params.box[2] = static_cast<float>(pme.box[2]);
                        stored_pme_params.alpha = pme.alpha;
                        stored_pme_params.set = true;
                        
                        py::print("Debug: Read from global PME parameters:",
                                  "cutoff=", pme.cutoff, 
                                  "box=[", pme.box[0], ",", pme.box[1], ",", pme.box[2], "]",
                                  "alpha=", pme.alpha);
                    }
                } catch (...) {
                    // Ignore any exceptions
                }
                
                if (!pme_initialized_globally) {
                    throw std::runtime_error("PME not properly initialized. Call initializePMEParameters first.");
                }
                
                if (!stored_pgp_params.set) {
                    throw std::runtime_error("PGP parameters not set. Call setPGPParameters first.");
                }
                
                // Create PGPContext
                try {
                    global_pgp_context = std::make_unique<::platform::cpu::energy::pgp::PGPContext>();
                } catch (const std::exception& e) {
                    throw std::runtime_error(std::string("Failed to create PGPContext: ") + e.what());
                }
                
                // Convert parameters
                std::array<double, 3> box_array = {
                    static_cast<double>(stored_pme_params.box[0]), 
                    static_cast<double>(stored_pme_params.box[1]), 
                    static_cast<double>(stored_pme_params.box[2])
                };
                std::array<int, 3> meshSize_array = {
                    stored_pgp_params.meshSize[0],
                    stored_pgp_params.meshSize[1],
                    stored_pgp_params.meshSize[2]
                };
                std::array<int, 3> gridSize_array = {
                    stored_pgp_params.potentialGridSize[0],
                    stored_pgp_params.potentialGridSize[1],
                    stored_pgp_params.potentialGridSize[2]
                };
                
                // Initialize PGPContext
                global_pgp_context->initialize(
                    stored_pme_params.cutoff,
                    box_array,
                    stored_pgp_params.alpha,
                    meshSize_array,
                    stored_pgp_params.potential_cutoff,
                    gridSize_array,
                    stored_pgp_params.splineOrder,
                    stored_pgp_params.tolerance
                );
                
                pgp_initialized = true;
            }
            
            // Use PGPContext for precomputation
            // Note: PGPContext's precompute_grid_potential works per atom type
            // We'll precompute for all atom types in the system
            for (int atom_type = 0; atom_type < state.forcefield.numTotalTypes; ++atom_type) {
                global_pgp_context->precomputeGridPotential(state, atom_type);
            }
            
            // Also call the original function to ensure global state is set for PGP Complete
            ::pygcmc::simulation::Simulation::precomputeGridPotential(state, fixed_only);
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
            std::lock_guard<std::mutex> lock(pgp_mutex);
            
            if (!pgp_initialized || !global_pgp_context) {
                // Fall back to original function
                double energy = 0.0;
                ::pygcmc::simulation::Simulation::interpolateMoleculeEnergy(state, energy);
                return energy;
            }
            
            // For now, fall back to original function since PGPContext doesn't have interpolateMoleculeEnergy
            double energy = 0.0;
            ::pygcmc::simulation::Simulation::interpolateMoleculeEnergy(state, energy);
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
            std::lock_guard<std::mutex> lock(pgp_mutex);
            
            if (!pgp_initialized || !global_pgp_context) {
                // Fall back to original function
                return ::pygcmc::simulation::Simulation::calculateMoleculeEnergy(state);
            }
            
            // For now, fall back to original function since PGPContext doesn't have interpolateMoleculeEnergy
            return ::pygcmc::simulation::Simulation::calculateMoleculeEnergy(state);
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
            std::lock_guard<std::mutex> lock(pgp_mutex);
            
            if (!pgp_initialized || !global_pgp_context) {
                throw std::runtime_error("PGP not initialized. Call initializePMEParameters first.");
            }
            
            // Use PGPContext
            auto energy = global_pgp_context->computeSystemEnergy(state);
            
            // Store results in state for compatibility
            state.ewald_energy.real_space = energy.real_space;
            state.ewald_energy.reciprocal = energy.reciprocal;
            state.ewald_energy.self = energy.self;
            state.ewald_energy.total = energy.total;
            
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
            std::lock_guard<std::mutex> lock(pgp_mutex);
            
            if (!pgp_initialized || !global_pgp_context) {
                throw std::runtime_error("PGP not initialized. Call initializePMEParameters first.");
            }
            
            // Use PGPContext for movement energy
                // Extract movement residue indices
                std::vector<int> movementIndices;
                for(const auto& movementInfo : state.movementResidues) {
                    for(int i = movementInfo.startIndex;
                        i < movementInfo.startIndex + movementInfo.activeCount; i++) {
                        movementIndices.push_back(i);
                    }
                }
                auto energy = global_pgp_context->computeMovementEnergy(state, movementIndices);
                
                // Store results in state for compatibility
                state.ewald_energy.real_space = energy.real_space;
                state.ewald_energy.reciprocal = energy.reciprocal;
                state.ewald_energy.self = energy.self;
                state.ewald_energy.total = energy.total;
            
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
            ::pygcmc::simulation::Simulation::computeSystemEnergyPGPComplete(state);
            
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
            ::pygcmc::simulation::Simulation::computeMovementEnergyPGPComplete(state);
            
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
}

} // namespace simulation
} // namespace bindings
} // namespace pygcmc