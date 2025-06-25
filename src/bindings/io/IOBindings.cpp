// src/bindings/io/IOBindings.cpp

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

void init_topology_bindings(py::module& m, py::module& io_module) {
    // PSFParser bindings
    auto psf_parser = py::class_<pygcmc::io::PSFParser>(m, "PSFParser")
        .def(py::init<>())
        .def("parse_to_topology", &pygcmc::io::PSFParser::parse_to_topology)
        .def_static("parse_file", &pygcmc::io::PSFParser::parse_file,
            py::arg("filename"),
            "Parse a PSF file and return a new Topology object")
        .def_static("parse_string", &pygcmc::io::PSFParser::parse_string,
            py::arg("psf_str"),
            "Parse a PSF string and return a new Topology object");
    io_module.attr("PSFParser") = psf_parser;

    // TOPParser bindings
    auto top_parser = py::class_<pygcmc::io::TOPParser>(m, "TOPParser")
        .def(py::init<>())
        .def("parse_to_topology", &pygcmc::io::TOPParser::parse_to_topology)
        .def_static("parse_file", &pygcmc::io::TOPParser::parse_file,
            py::arg("filename"),
            "Parse a topology file and return a new Topology object")
        .def_static("parse_string", &pygcmc::io::TOPParser::parse_string,
            py::arg("top_str"),
            "Parse a topology string and return a new Topology object")
        .def_static("enable_debug", &pygcmc::io::TOPParser::enable_debug,
            py::arg("enable"),
            "Enable or disable debug output")
        .def_static("is_debug_enabled", &pygcmc::io::TOPParser::is_debug_enabled,
            "Check if debug output is enabled");
    io_module.attr("TOPParser") = top_parser;
}

void init_forcefield_bindings(py::module& m, py::module& io_module) {
    // PRMParser bindings
    auto prm_parser = py::class_<pygcmc::io::PRMParser>(m, "PRMParser")
        .def(py::init<>())
        .def_static("parse_string", &pygcmc::io::PRMParser::parse_string,
            py::arg("content"), py::arg("ff"))
        .def_static("parse_file", &pygcmc::io::PRMParser::parse_file,
            py::arg("filename"))
        .def_static("parse_files", &pygcmc::io::PRMParser::parse_files,
            py::arg("filenames"))
        .def_static("parse_file_to_forcefield", &pygcmc::io::PRMParser::parse_file_to_forcefield,
            py::arg("filename"), py::arg("ff"))
        .def("parse", &pygcmc::io::PRMParser::parse,
            py::arg("filename"), py::arg("ff"));

    // Add to both main module and io submodule for backward compatibility
    m.attr("PRMParser") = prm_parser;
    io_module.attr("PRMParser") = prm_parser;
}

void init_parameters_bindings(py::module& m, py::module& io_module) {
    // INPParser bindings
    auto inp_parser = py::class_<pygcmc::io::INPParser>(io_module, "INPParser")
        .def_static("parse_file", &pygcmc::io::INPParser::parse_file,
            py::arg("filename"),
            "Parse an input file and return a new Param object")
        .def_static("parse_string", &pygcmc::io::INPParser::parse_string,
            py::arg("content"),
            "Parse an input string and return a new Param object")
        .def_static("parse_to_param", &pygcmc::io::INPParser::parse_to_param,
            py::arg("filename"), py::arg("param"),
            "Parse an input file into an existing Param object")
        .def_static("parse_string_to_param", &pygcmc::io::INPParser::parse_string_to_param,
            py::arg("content"), py::arg("param"),
            "Parse an input string into an existing Param object");
    
    // Add to main module for backward compatibility
    m.attr("INPParser") = inp_parser;
}

void init_io_bindings(py::module& m) {
    // Create io submodule
    auto io_module = m.def_submodule("io", "Input/Output operations");
    
    // Initialize all IO binding groups
    init_structure_bindings(m, io_module);
    init_topology_bindings(m, io_module);
    init_forcefield_bindings(m, io_module);
    init_parameters_bindings(m, io_module);
}

} // namespace io
} // namespace bindings
} // namespace pygcmc