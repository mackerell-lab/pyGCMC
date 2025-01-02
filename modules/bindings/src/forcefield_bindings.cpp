#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pygcmc/core/forcefield.hpp"
#include "pygcmc/core/io/ff_parser.hpp"

namespace py = pybind11;
using namespace pygcmc::core;

void init_forcefield_bindings(py::module& m) {
    // Bind ForceField class
    py::class_<ForceField, std::shared_ptr<ForceField>>(m, "ForceField")
        .def(py::init<>())
        .def_property("cutoff", &ForceField::get_cutoff, &ForceField::set_cutoff)
        .def_property("switching", &ForceField::get_switching, &ForceField::set_switching)
        .def_property("pairlist_distance", &ForceField::get_pairlist_distance, &ForceField::set_pairlist_distance)
        .def_property_readonly("nonbonded_params", 
            static_cast<const std::map<std::string, io::ForceFieldPair>& (ForceField::*)() const>(&ForceField::nonbonded_params))
        .def_property_readonly("nbfix_params", 
            static_cast<const std::map<std::pair<std::string, std::string>, io::ForceFieldPair>& (ForceField::*)() const>(&ForceField::nbfix_params))
        .def("print_nonbonded_params", &ForceField::print_nonbonded_params, "Print all nonbonded parameters")
        // Add new data access methods
        .def("get_nonbonded_parameters", [](ForceField& self) {
            py::list result;
            for (const auto& [type, params] : self.nonbonded_params()) {
                py::dict param_data;
                param_data["atom_type"] = type;
                param_data["epsilon"] = params.epsilon;
                param_data["rmin"] = params.rmin;
                param_data["epsilon14"] = params.epsilon * 0.5;  // Example scaling
                param_data["rmin14"] = params.rmin;
                result.append(param_data);
            }
            return result;
        })
        .def("get_nbfix_parameters", [](ForceField& self) {
            py::list result;
            for (const auto& [types, params] : self.nbfix_params()) {
                py::dict param_data;
                param_data["type1"] = types.first;
                param_data["type2"] = types.second;
                param_data["epsilon"] = params.epsilon;
                param_data["rmin"] = params.rmin;
                result.append(param_data);
            }
            return result;
        })
        .def("get_global_parameters", [](ForceField& self) {
            py::dict params;
            params["cutoff"] = self.get_cutoff();
            params["switching"] = self.get_switching();
            params["pairlist_distance"] = self.get_pairlist_distance();
            return params;
        });
} 