// src/bindings/system/SystemMolecularBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "system/SystemModule.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace system {

void init_molecular_bindings(py::module& m, py::module&) {
    py::class_<::pygcmc::system::MolecularSystem, std::shared_ptr<::pygcmc::system::MolecularSystem>>(m, "MolecularSystem")
        .def(py::init<>())
        .def("combine", 
             [](::pygcmc::system::MolecularSystem& self, 
                py::object structure, 
                py::object topology) {
                 // First check for None values
                 if (structure.is_none() || topology.is_none()) {
                     throw py::value_error("Structure and Topology cannot be None");
                 }
                 
                 // Then try to convert to the correct types
                 const pygcmc::model::Structure* struct_ptr = structure.cast<const pygcmc::model::Structure*>();
                 const pygcmc::model::Topology* top_ptr = topology.cast<const pygcmc::model::Topology*>();
                 
                 // Create shared_ptr without ownership
                 auto struct_shared = std::shared_ptr<pygcmc::model::Structure>(
                     const_cast<pygcmc::model::Structure*>(struct_ptr), 
                     [](pygcmc::model::Structure*){});
                 auto top_shared = std::shared_ptr<pygcmc::model::Topology>(
                     const_cast<pygcmc::model::Topology*>(top_ptr), 
                     [](pygcmc::model::Topology*){});
                 
                 return self.combine(struct_shared, top_shared);
             },
             py::arg("structure"), 
             py::arg("topology"),
             py::return_value_policy::move,
             "Combine Structure and Topology data into a Molecular object")
        .def("combine_multiple",
             [](::pygcmc::system::MolecularSystem& self,
                py::object structure,
                py::list topologies) {
                 // Check for None values
                 if (structure.is_none() || topologies.is_none()) {
                     throw py::value_error("Structure and topologies list cannot be None");
                 }
                 
                 // Convert structure to shared_ptr
                 const pygcmc::model::Structure* struct_ptr = structure.cast<const pygcmc::model::Structure*>();
                 auto struct_shared = std::shared_ptr<pygcmc::model::Structure>(
                     const_cast<pygcmc::model::Structure*>(struct_ptr),
                     [](pygcmc::model::Structure*){});
                 
                 // Convert list of topologies to vector of shared_ptr
                 std::vector<std::shared_ptr<pygcmc::model::Topology>> top_vec;
                 for (const auto& top : topologies) {
                     const pygcmc::model::Topology* top_ptr = top.cast<const pygcmc::model::Topology*>();
                     top_vec.push_back(std::shared_ptr<pygcmc::model::Topology>(
                         const_cast<pygcmc::model::Topology*>(top_ptr),
                         [](pygcmc::model::Topology*){}));
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
              const pygcmc::model::Structure* struct_ptr = structure.cast<const pygcmc::model::Structure*>();
              const pygcmc::model::Topology* top_ptr = topology.cast<const pygcmc::model::Topology*>();
              
              // Create shared_ptr without ownership
              auto struct_shared = std::shared_ptr<pygcmc::model::Structure>(
                  const_cast<pygcmc::model::Structure*>(struct_ptr), 
                  [](pygcmc::model::Structure*){});
              auto top_shared = std::shared_ptr<pygcmc::model::Topology>(
                  const_cast<pygcmc::model::Topology*>(top_ptr), 
                  [](pygcmc::model::Topology*){});
              
              // Create MolecularSystem and combine
              auto mol_system = std::make_shared<::pygcmc::system::MolecularSystem>();
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
              const pygcmc::model::Structure* struct_ptr = structure.cast<const pygcmc::model::Structure*>();
              auto struct_shared = std::shared_ptr<pygcmc::model::Structure>(
                  const_cast<pygcmc::model::Structure*>(struct_ptr), 
                  [](pygcmc::model::Structure*){});
              
              // Convert topologies to vector of shared_ptr
              std::vector<std::shared_ptr<pygcmc::model::Topology>> top_vec;
              for (const auto& top : topologies) {
                  const pygcmc::model::Topology* top_ptr = top.cast<const pygcmc::model::Topology*>();
                  top_vec.push_back(std::shared_ptr<pygcmc::model::Topology>(
                      const_cast<pygcmc::model::Topology*>(top_ptr),
                      [](pygcmc::model::Topology*){}));
              }
              
              // Create MolecularSystem and combine
              auto mol_system = std::make_shared<::pygcmc::system::MolecularSystem>();
              return mol_system->combine_multiple(struct_shared, top_vec);
          },
          py::return_value_policy::move,
          "Directly combine Structure with multiple Topology files into a Molecular object");
}

} // namespace system
} // namespace bindings
} // namespace pygcmc