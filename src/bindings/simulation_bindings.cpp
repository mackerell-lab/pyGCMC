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
}

} // namespace bindings
} // namespace pygcmc
