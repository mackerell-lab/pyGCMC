// src/bindings/system_bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "system/SystemModule.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {

void init_system(py::module& m) {
    // Add LogLevel enum
    py::enum_<system::common::LogLevel>(m, "LogLevel")
        .value("DEBUG", system::common::LogLevel::DEBUG)
        .value("INFO", system::common::LogLevel::INFO)
        .value("WARNING", system::common::LogLevel::WARNING)
        .value("ERROR", system::common::LogLevel::ERROR);

    py::class_<system::System>(m, "System")
        .def(py::init<>())
        .def_static("set_verbose", &system::System::set_verbose, 
                   "Set verbose mode for logging")
        .def_static("set_log_level", &system::System::set_log_level, 
                   "Set the minimum log level")
        .def("initialize_parameters", &system::System::initialize_parameters, 
             "Initialize all system parameters")
        .def("process_cavity_list", &system::System::process_cavity_list,
             "Process and sort the cavity list")
        .def("initialize_mc_time_list", &system::System::initialize_mc_time_list,
             "Initialize Monte Carlo time list");

    py::class_<system::MolecularSystem, std::shared_ptr<system::MolecularSystem>>(m, "MolecularSystem")
        .def(py::init<>())
        .def("combine", 
             [](system::MolecularSystem& self, 
                py::object structure, 
                py::object topology) {
                 // First check for None values
                 if (structure.is_none() || topology.is_none()) {
                     throw py::value_error("Structure and Topology cannot be None");
                 }
                 
                 // Then try to convert to the correct types
                 const model::Structure* struct_ptr = structure.cast<const model::Structure*>();
                 const model::Topology* top_ptr = topology.cast<const model::Topology*>();
                 
                 // Create shared_ptr without ownership
                 auto struct_shared = std::shared_ptr<model::Structure>(
                     const_cast<model::Structure*>(struct_ptr), 
                     [](model::Structure*){});
                 auto top_shared = std::shared_ptr<model::Topology>(
                     const_cast<model::Topology*>(top_ptr), 
                     [](model::Topology*){});
                 
                 return self.combine(struct_shared, top_shared);
             },
             py::arg("structure"), 
             py::arg("topology"),
             py::return_value_policy::move,
             "Combine Structure and Topology data into a Molecular object")
        .def("combine_multiple",
             [](system::MolecularSystem& self,
                py::object structure,
                py::list topologies) {
                 // Check for None values
                 if (structure.is_none() || topologies.is_none()) {
                     throw py::value_error("Structure and topologies list cannot be None");
                 }
                 
                 // Convert structure to shared_ptr
                 const model::Structure* struct_ptr = structure.cast<const model::Structure*>();
                 auto struct_shared = std::shared_ptr<model::Structure>(
                     const_cast<model::Structure*>(struct_ptr),
                     [](model::Structure*){});
                 
                 // Convert list of topologies to vector of shared_ptr
                 std::vector<std::shared_ptr<model::Topology>> top_vec;
                 for (const auto& top : topologies) {
                     const model::Topology* top_ptr = top.cast<const model::Topology*>();
                     top_vec.push_back(std::shared_ptr<model::Topology>(
                         const_cast<model::Topology*>(top_ptr),
                         [](model::Topology*){}));
                 }
                 
                 return self.combine_multiple(struct_shared, top_vec);
             },
             py::arg("structure"),
             py::arg("topologies"),
             py::return_value_policy::move,
             "Combine Structure with multiple Topology files into a Molecular object");

    // Add direct Combine function
    m.def("Combine", 
          [](py::object structure, 
             py::object topology) {
              // First check for None values
              if (structure.is_none() || topology.is_none()) {
                  throw py::value_error("Structure and Topology cannot be None");
              }
              
              // Then try to convert to the correct types
              const model::Structure* struct_ptr = structure.cast<const model::Structure*>();
              const model::Topology* top_ptr = topology.cast<const model::Topology*>();
              
              // Create shared_ptr without ownership
              auto struct_shared = std::shared_ptr<model::Structure>(
                  const_cast<model::Structure*>(struct_ptr), 
                  [](model::Structure*){});
              auto top_shared = std::shared_ptr<model::Topology>(
                  const_cast<model::Topology*>(top_ptr), 
                  [](model::Topology*){});
              
              // Create MolecularSystem and combine
              auto mol_system = std::make_shared<system::MolecularSystem>();
              return mol_system->combine(struct_shared, top_shared);
          },
          py::arg("structure"), 
          py::arg("topology"),
          py::return_value_policy::move,
          "Directly combine Structure and Topology data into a Molecular object");

    // Add direct Combine function with multiple topologies
    m.def("Combine", 
          [](py::object structure, 
             py::args topologies) {
              // First check for None values
              if (structure.is_none() || topologies.empty()) {
                  throw py::value_error("Structure and at least one Topology must be provided");
              }
              
              // Convert structure to shared_ptr
              const model::Structure* struct_ptr = structure.cast<const model::Structure*>();
              auto struct_shared = std::shared_ptr<model::Structure>(
                  const_cast<model::Structure*>(struct_ptr), 
                  [](model::Structure*){});
              
              // Convert topologies to vector of shared_ptr
              std::vector<std::shared_ptr<model::Topology>> top_vec;
              for (const auto& top : topologies) {
                  const model::Topology* top_ptr = top.cast<const model::Topology*>();
                  top_vec.push_back(std::shared_ptr<model::Topology>(
                      const_cast<model::Topology*>(top_ptr),
                      [](model::Topology*){}));
              }
              
              // Create MolecularSystem and combine
              auto mol_system = std::make_shared<system::MolecularSystem>();
              return mol_system->combine_multiple(struct_shared, top_vec);
          },
          py::return_value_policy::move,
          "Directly combine Structure with multiple Topology files into a Molecular object");

    // Bind TypeMaps
    py::class_<pygcmc::model::TypeMaps>(m, "TypeMaps")
        .def(py::init<>())
        .def("get_or_add_type", &pygcmc::model::TypeMaps::getOrAddType)
        .def("get_type_name", &pygcmc::model::TypeMaps::getTypeName)
        .def_readonly("atomTypes", &pygcmc::model::TypeMaps::atomTypes);

    // Bind MonteCarloSystem
    py::class_<pygcmc::system::MonteCarloSystem>(m, "MonteCarloSystem")
        .def(py::init<>())
        .def("initialize", &pygcmc::system::MonteCarloSystem::initialize)
        .def("set_force_field", &pygcmc::system::MonteCarloSystem::setForceField)
        .def("initialize_force_field", &pygcmc::system::MonteCarloSystem::initializeForceField)
        .def("initialize_from_molecular", [](pygcmc::system::MonteCarloSystem& self, py::object molecular) {
            if (molecular.is_none()) {
                throw py::value_error("Molecular object cannot be None");
            }
            
            try {
                // First try MolecularSystem
                auto* molSys = molecular.cast<pygcmc::system::MolecularSystem*>();
                if (molSys) {
                    self.initializeFromMolecular(molSys->get_molecular());
                    return;
                }
            } catch (py::cast_error&) {}
            
            try {
                // Then try Molecular directly
                auto mol = molecular.cast<std::shared_ptr<pygcmc::model::Molecular>>();
                if (mol) {
                    self.initializeFromMolecular(mol);
                    return;
                }
            } catch (py::cast_error&) {}
            
            throw py::type_error("Argument must be either MolecularSystem or Molecular");
        })
        .def("add_movement_molecules", [](pygcmc::system::MonteCarloSystem& self, py::list molecules) {
            std::vector<pygcmc::system::montecarlo::MovementMolecularInfo> mol_vec;
            for (const auto& mol : molecules) {
                mol_vec.push_back(mol.cast<pygcmc::system::montecarlo::MovementMolecularInfo>());
            }
            self.addMovementMolecules(mol_vec);
        }, py::arg("molecules"), "Add movement molecules for GCMC simulation")
        .def("get_type_maps", &pygcmc::system::MonteCarloSystem::getTypeMaps, py::return_value_policy::reference)
        .def("insert_residue", &pygcmc::system::MonteCarloSystem::insertResidue)
        .def("remove_residue", &pygcmc::system::MonteCarloSystem::removeResidue)
        .def("translate_residue", &pygcmc::system::MonteCarloSystem::translateResidue)
        .def("calc_non_bonded_energy", &pygcmc::system::MonteCarloSystem::calcNonBondedEnergy)
        .def("calc_total_energy", &pygcmc::system::MonteCarloSystem::calcTotalEnergy)
        .def("get_state", (const pygcmc::model::MCState& (pygcmc::system::MonteCarloSystem::*)() const) &pygcmc::system::MonteCarloSystem::getState, py::return_value_policy::reference)
        .def("get_state_mutable", (pygcmc::model::MCState& (pygcmc::system::MonteCarloSystem::*)()) &pygcmc::system::MonteCarloSystem::getState, py::return_value_policy::reference)
        .def("get_active_atom_count", &pygcmc::system::MonteCarloSystem::getActiveAtomCount)
        .def("get_active_residue_count", &pygcmc::system::MonteCarloSystem::getActiveResidueCount)
        .def("set_switching_function", &pygcmc::system::MonteCarloSystem::setSwitchingFunction,
             py::arg("enable"), py::arg("r_on") = 1.0f, py::arg("r_off") = 1.2f,
             "Set or disable CHARMM-style smooth switching function")
        .def("calculate_switching_function", &pygcmc::system::MonteCarloSystem::calculateSwitchingFunction,
             py::arg("r"), "Calculate switching function value at the given distance")
        .def("is_using_switching_function", &pygcmc::system::MonteCarloSystem::isUsingSwitchingFunction,
             "Get whether the switching function is currently enabled")
        .def("get_switching_r_on", &pygcmc::system::MonteCarloSystem::getSwitchingROn,
             "Get the inner cutoff radius")
        .def("get_switching_r_off", &pygcmc::system::MonteCarloSystem::getSwitchingROff,
             "Get the outer cutoff radius")
        .def("apply_switching_to_state", &pygcmc::system::MonteCarloSystem::applySwitchingToState,
             py::arg("state"), "Apply the current switching function settings to an external state object");

    // Bind MovementMolecularInfo
    py::class_<pygcmc::system::montecarlo::MovementMolecularInfo>(m, "MovementMolecularInfo")
        .def(py::init<std::shared_ptr<pygcmc::model::Molecular>, int>(),
             py::arg("molecular"),
             py::arg("maxCopies"))
        .def_readwrite("molecular", &pygcmc::system::montecarlo::MovementMolecularInfo::molecular)
        .def_readwrite("maxCopies", &pygcmc::system::montecarlo::MovementMolecularInfo::maxCopies);
}

} // namespace bindings
} // namespace pygcmc
