// src/bindings/io_bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "io/pdbParser.hpp"
#include "io/psfParser.hpp"
#include "io/topParser.hpp"
#include "io/prmParser.hpp"
#include "io/inpParser.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {

void init_io(py::module& m) {
    // Create io submodule
    auto io = m.def_submodule("io", "Input/Output operations");
    
    // Bind PDBParser to io submodule
    auto pdb_parser = py::class_<io::PDBParser>(io, "PDBParser")
        .def_static("parse_file", &io::PDBParser::parse_file,
            py::arg("filename"),
            "Parse PDB file and return Structure object")
        .def_static("parse_string", &io::PDBParser::parse_string,
            py::arg("pdbStr"),
            "Parse PDB string and return Structure object")
        .def_static("parse_to_structure", &io::PDBParser::parse_to_structure,
            py::arg("filename"), py::arg("structure"),
            "Parse PDB file and populate Structure object")
        .def_static("parse_string_to_structure", &io::PDBParser::parse_string_to_structure,
            py::arg("pdbStr"), py::arg("structure"),
            "Parse PDB string and populate Structure object");
    m.attr("PDBParser") = pdb_parser;  // 将 PDBParser 也添加到主模块

    // PSFParser bindings
    auto psf_parser = py::class_<io::PSFParser>(m, "PSFParser")
        .def(py::init<>())
        .def("parse_to_topology", &io::PSFParser::parse_to_topology)
        .def_static("parse_file", &io::PSFParser::parse_file,
            py::arg("filename"),
            "Parse a PSF file and return a new Topology object")
        .def_static("parse_string", &io::PSFParser::parse_string,
            py::arg("psf_str"),
            "Parse a PSF string and return a new Topology object");
    io.attr("PSFParser") = psf_parser;

    // TOPParser bindings
    auto top_parser = py::class_<io::TOPParser>(m, "TOPParser")
        .def(py::init<>())
        .def("parse_to_topology", &io::TOPParser::parse_to_topology)
        .def_static("parse_file", &io::TOPParser::parse_file,
            py::arg("filename"),
            "Parse a topology file and return a new Topology object")
        .def_static("parse_string", &io::TOPParser::parse_string,
            py::arg("top_str"),
            "Parse a topology string and return a new Topology object")
        .def_static("enable_debug", &io::TOPParser::enable_debug,
            py::arg("enable"),
            "Enable or disable debug output")
        .def_static("is_debug_enabled", &io::TOPParser::is_debug_enabled,
            "Check if debug output is enabled");
    io.attr("TOPParser") = top_parser;

    // Add PRMParser bindings
    auto prm_parser = py::class_<PRMParser>(m, "PRMParser")
        // Static methods
        .def_static("parse_string", &PRMParser::parse_string,
            "Parse parameter content from a string",
            py::arg("content"), py::arg("ff"))
        .def_static("parse_file", &PRMParser::parse_file,
            "Parse parameter content from a file and return a new ForceField object",
            py::arg("filename"))
        .def_static("parse_files", &PRMParser::parse_files,
            "Parse parameter content from multiple files and return a new ForceField object",
            py::arg("filenames"))
        .def_static("parse_file_to_forcefield", &PRMParser::parse_file_to_forcefield,
            "Parse parameter content from a file into an existing ForceField object",
            py::arg("filename"), py::arg("ff"))
        // Instance methods
        .def("parse", &PRMParser::parse,
            "Parse parameter content from a file (instance method)",
            py::arg("filename"), py::arg("ff"));

    m.attr("PRMParser") = prm_parser;  // Add to main module as well
    io.attr("PRMParser") = prm_parser;  // Add to io submodule

    // Add InpParser bindings
    auto inp_parser = py::class_<io::InpParser>(io, "InpParser")
        .def_static("parse_file", &io::InpParser::parse_file,
            py::arg("filename"),
            "Parse an input file and return a new Param object")
        .def_static("parse_string", &io::InpParser::parse_string,
            py::arg("content"),
            "Parse an input string and return a new Param object")
        .def_static("parse_to_param", &io::InpParser::parse_to_param,
            py::arg("filename"), py::arg("param"),
            "Parse an input file into an existing Param object")
        .def_static("parse_string_to_param", &io::InpParser::parse_string_to_param,
            py::arg("content"), py::arg("param"),
            "Parse an input string into an existing Param object");
    m.attr("InpParser") = inp_parser;  // Add to main module as well
}

} // namespace bindings
} // namespace pygcmc


