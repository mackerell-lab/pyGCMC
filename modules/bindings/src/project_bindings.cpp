// modules/bindings/src/project_bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pygcmc/core/project.hpp"

namespace py = pybind11;
using namespace pygcmc::core;

void init_project_bindings(py::module& m) {
    // Bind Project class
    py::class_<Project>(m, "Project")
        .def(py::init<const std::string&>(), py::arg("name") = "")
        .def("create_structure", &Project::create_structure, "Create a new empty structure")
        .def("load_structure", &Project::load_structure, 
             py::arg("pdb_file"), 
             py::arg("top_file") = "",
             "Load structure from PDB and topology files")
        .def("load_forcefield", &Project::load_forcefield, "Load force field from parameter files")
        .def("get_name", &Project::get_name, "Get project name")
        .def("print_atom_info", &Project::print_atom_info, "Print detailed information for a single atom")
        .def("print_detailed_atom_info", &Project::print_detailed_atom_info, "Print detailed information for a single atom in table format")
        .def("print_atom_table_header", &Project::print_atom_table_header, "Print header for atom table")
        .def("print_all_atoms", &Project::print_all_atoms, "Print information for all atoms")
        .def("print_forcefield_info", &Project::print_forcefield_info, "Print force field information")
        .def("print_nbfix_info", &Project::print_nbfix_info, "Print NBFIX parameters")
        .def("print_global_parameters", &Project::print_global_parameters, "Print global force field parameters");
} 