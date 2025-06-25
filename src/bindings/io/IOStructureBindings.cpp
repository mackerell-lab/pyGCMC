// src/bindings/io/IOStructureBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "io/IOModule.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace io {

void init_structure_bindings(py::module& m, py::module& io_module) {
    // PDBParser bindings for structure parsing
    auto pdb_parser = py::class_<pygcmc::io::PDBParser>(io_module, "PDBParser")
        .def_static("parse_file", &pygcmc::io::PDBParser::parse_file,
            py::arg("filename"),
            "Parse PDB file and return Structure object")
        .def_static("parse_string", &pygcmc::io::PDBParser::parse_string,
            py::arg("pdbStr"),
            "Parse PDB string and return Structure object")
        .def_static("parse_to_structure", &pygcmc::io::PDBParser::parse_to_structure,
            py::arg("filename"), py::arg("structure"),
            "Parse PDB file and populate Structure object")
        .def_static("parse_string_to_structure", &pygcmc::io::PDBParser::parse_string_to_structure,
            py::arg("pdbStr"), py::arg("structure"),
            "Parse PDB string and populate Structure object");
    
    // Also add to main module for backward compatibility
    m.attr("PDBParser") = pdb_parser;
}

} // namespace io
} // namespace bindings
} // namespace pygcmc