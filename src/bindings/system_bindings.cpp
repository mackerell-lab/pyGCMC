#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "system/system.hpp"
#include "system/molecularSystem.hpp"
#include "system/MonteCarloSystem.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {

void init_system(py::module& m) {
    // Add LogLevel enum
    py::enum_<system::LogLevel>(m, "LogLevel")
        .value("DEBUG", system::LogLevel::DEBUG)
        .value("INFO", system::LogLevel::INFO)
        .value("WARNING", system::LogLevel::WARNING)
        .value("ERROR", system::LogLevel::ERROR);

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
    py::class_<gcmc::MonteCarloSystem::TypeMaps>(m, "TypeMaps")
        .def("get_or_add_type", &gcmc::MonteCarloSystem::TypeMaps::getOrAddType)
        .def("get_type_name", &gcmc::MonteCarloSystem::TypeMaps::getTypeName)
        .def_readonly("atom_types", &gcmc::MonteCarloSystem::TypeMaps::atomTypes);

    // Bind MonteCarloSystem
    py::class_<gcmc::MonteCarloSystem>(m, "MonteCarloSystem")
        .def(py::init<>())
        .def("initialize", &gcmc::MonteCarloSystem::initialize)
        .def("set_force_field", &gcmc::MonteCarloSystem::setForceField)
        .def("initialize_from_molecular", &gcmc::MonteCarloSystem::initializeFromMolecular)
        .def("get_type_maps", &gcmc::MonteCarloSystem::getTypeMaps, py::return_value_policy::reference)
        .def("insert_residue", &gcmc::MonteCarloSystem::insertResidue)
        .def("remove_residue", &gcmc::MonteCarloSystem::removeResidue)
        .def("translate_residue", &gcmc::MonteCarloSystem::translateResidue)
        .def("calc_non_bonded_energy", &gcmc::MonteCarloSystem::calcNonBondedEnergy)
        .def("calc_total_energy", &gcmc::MonteCarloSystem::calcTotalEnergy)
        .def("get_state", (const gcmc::SystemState& (gcmc::MonteCarloSystem::*)() const) &gcmc::MonteCarloSystem::getState, py::return_value_policy::reference)
        .def("get_active_atom_count", &gcmc::MonteCarloSystem::getActiveAtomCount)
        .def("get_active_residue_count", &gcmc::MonteCarloSystem::getActiveResidueCount);
}

} // namespace bindings
} // namespace pygcmc
