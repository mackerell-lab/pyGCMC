// src/bindings/io/IOStructureBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "io/IOModule.hpp"
#include "../../io/common/IOConfig.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace io {

void init_structure_bindings(py::module& m, py::module& io_module) {
    // Add IO configuration functions
    io_module.def("set_io_verbose", [](bool verbose) {
            pygcmc::io::IOConfig::verbose_errors = verbose;
        },
        py::arg("verbose") = false,
        "Set whether IO operations should print error messages (default: False)");

    io_module.def("get_io_verbose", []() {
            return pygcmc::io::IOConfig::verbose_errors;
        },
        "Get current IO verbosity setting");

    // IOConfig class (optional, for direct access)
    py::class_<pygcmc::io::IOConfig>(io_module, "IOConfig", "Configuration for IO operations")
        .def_property_static("verbose_errors",
            [](py::object) { return pygcmc::io::IOConfig::verbose_errors; },
            [](py::object, bool value) { pygcmc::io::IOConfig::verbose_errors = value; },
            "Control whether to print error messages to stderr (default: False)");

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
