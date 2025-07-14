// src/bindings/simulation/SimulationPGPIndependentBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "../../simulation/simulation.hpp"
#include "../../model/ModelModule.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace simulation {

void init_pgp_independent_bindings(py::module& m) {
    // Independent PGP parameter setting
    m.def("setPGPParametersIndependent",
        [](float alpha, const std::vector<int>& meshSize, float potential_cutoff,
           const std::vector<int>& potentialGridSize, int splineOrder, float tolerance) {
            if (meshSize.size() != 3 || potentialGridSize.size() != 3) {
                throw std::runtime_error("meshSize and potentialGridSize must have exactly three elements");
            }
            
            int meshSize_array[3] = {meshSize[0], meshSize[1], meshSize[2]};
            int potentialGridSize_array[3] = {potentialGridSize[0], potentialGridSize[1], potentialGridSize[2]};
            
            ::pygcmc::simulation::Simulation::setPGPParametersIndependent(
                static_cast<double>(alpha), meshSize_array, static_cast<double>(potential_cutoff), 
                potentialGridSize_array, splineOrder, static_cast<double>(tolerance));
        },
        "Set PGP parameters independently (no PME coupling)",
        py::arg("alpha"),
        py::arg("meshSize"),
        py::arg("potential_cutoff"),
        py::arg("potentialGridSize"),
        py::arg("splineOrder") = 4,
        py::arg("tolerance") = 1e-5f);
    
    // Independent PGP initialization
    m.def("initializePGPParametersIndependent",
        [](float cutoff, const std::vector<float>& box, float alpha,
           const std::vector<int>& meshSize, float potential_cutoff,
           const std::vector<int>& potentialGridSize, int splineOrder, float tolerance) {
            if (box.size() != 3) {
                throw std::runtime_error("box must have exactly three elements");
            }
            
            // Use double array to match C++ function signature
            double box_array[3] = {static_cast<double>(box[0]), 
                                   static_cast<double>(box[1]), 
                                   static_cast<double>(box[2])};
            
            // Handle optional parameters
            if (meshSize.empty() || potentialGridSize.empty()) {
                ::pygcmc::simulation::Simulation::initializePGPParametersIndependent(
                    static_cast<double>(cutoff), box_array, static_cast<double>(alpha), 
                    nullptr, static_cast<double>(potential_cutoff), 
                    nullptr, splineOrder, static_cast<double>(tolerance));
            } else {
                if (meshSize.size() != 3 || potentialGridSize.size() != 3) {
                    throw std::runtime_error("meshSize and potentialGridSize must have exactly three elements");
                }
                int meshSize_array[3] = {meshSize[0], meshSize[1], meshSize[2]};
                int potentialGridSize_array[3] = {potentialGridSize[0], potentialGridSize[1], potentialGridSize[2]};
                
                ::pygcmc::simulation::Simulation::initializePGPParametersIndependent(
                    static_cast<double>(cutoff), box_array, static_cast<double>(alpha), 
                    meshSize_array, static_cast<double>(potential_cutoff),
                    potentialGridSize_array, splineOrder, static_cast<double>(tolerance));
            }
        },
        "Initialize PGP parameters independently with automatic optimization",
        py::arg("cutoff"),
        py::arg("box"),
        py::arg("alpha") = 0.0f,
        py::arg("meshSize") = std::vector<int>(),
        py::arg("potential_cutoff") = 0.0f,
        py::arg("potentialGridSize") = std::vector<int>(),
        py::arg("splineOrder") = 4,
        py::arg("tolerance") = 1e-5f);
    
    // Independent energy calculation functions
    m.def("computeSystemEnergyPGPIndependent",
        [](::pygcmc::model::MCState& state) {
            ::pygcmc::simulation::Simulation::computeSystemEnergyPGPIndependent(state);
            
            // Return energy components
            py::dict energy_dict;
            energy_dict["real_space"] = state.ewald_energy.real_space;
            energy_dict["reciprocal"] = state.ewald_energy.reciprocal;
            energy_dict["self"] = state.ewald_energy.self;
            energy_dict["total"] = state.ewald_energy.total;
            
            return energy_dict;
        },
        "Calculate system energy using independent PGP implementation");
    
    m.def("computeMovementEnergyPGPIndependent",
        [](::pygcmc::model::MCState& state) {
            ::pygcmc::simulation::Simulation::computeMovementEnergyPGPIndependent(state);
            
            // Return energy components
            py::dict energy_dict;
            energy_dict["real_space"] = state.ewald_energy.real_space;
            energy_dict["reciprocal"] = state.ewald_energy.reciprocal;
            energy_dict["self"] = state.ewald_energy.self;
            energy_dict["total"] = state.ewald_energy.total;
            
            return energy_dict;
        },
        "Calculate movement energy using independent PGP implementation");
}

} // namespace simulation
} // namespace bindings
} // namespace pygcmc