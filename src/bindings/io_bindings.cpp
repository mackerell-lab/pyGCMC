// src/bindings/io_bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../io/pdbParser.hpp"
#include "../io/psfParser.hpp"
#include "../io/topParser.hpp"
#include "model/atom.hpp"
#include "model/residue.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {

void init_io(py::module& m) {
    // Create io submodule
    auto io = m.def_submodule("io", "Input/Output operations");
    
    // Bind HelixInfo struct in both main module and io submodule
    auto helix_info = py::class_<pygcmc::io::PDBParser::HelixInfo>(m, "HelixInfo")
        .def(py::init<>())
        .def_readwrite("helixId", &pygcmc::io::PDBParser::HelixInfo::helixId)
        .def_readwrite("initResName", &pygcmc::io::PDBParser::HelixInfo::initResName)
        .def_readwrite("initChainId", &pygcmc::io::PDBParser::HelixInfo::initChainId)
        .def_readwrite("initSeqNum", &pygcmc::io::PDBParser::HelixInfo::initSeqNum)
        .def_readwrite("initICode", &pygcmc::io::PDBParser::HelixInfo::initICode)
        .def_readwrite("endResName", &pygcmc::io::PDBParser::HelixInfo::endResName)
        .def_readwrite("endChainId", &pygcmc::io::PDBParser::HelixInfo::endChainId)
        .def_readwrite("endSeqNum", &pygcmc::io::PDBParser::HelixInfo::endSeqNum)
        .def_readwrite("endICode", &pygcmc::io::PDBParser::HelixInfo::endICode)
        .def_readwrite("helixClass", &pygcmc::io::PDBParser::HelixInfo::helixClass);
    io.attr("HelixInfo") = helix_info;

    // Bind ParseResult struct in both main module and io submodule
    auto parse_result = py::class_<pygcmc::io::PDBParser::ParseResult>(m, "PDBParseResult")
        .def(py::init<>())
        .def_readwrite("atoms", &pygcmc::io::PDBParser::ParseResult::atoms)
        .def_readwrite("residues", &pygcmc::io::PDBParser::ParseResult::residues)
        .def_readwrite("terminals", &pygcmc::io::PDBParser::ParseResult::terminals)
        .def_readwrite("helices", &pygcmc::io::PDBParser::ParseResult::helices)
        .def_readwrite("sheets", &pygcmc::io::PDBParser::ParseResult::sheets)
        .def_readwrite("ssbonds", &pygcmc::io::PDBParser::ParseResult::ssbonds)
        .def_readwrite("boxDimensions", &pygcmc::io::PDBParser::ParseResult::boxDimensions);
    io.attr("PDBParseResult") = parse_result;

    // Bind PDBParser class in both main module and io submodule
    auto parser = py::class_<pygcmc::io::PDBParser>(m, "PDBParser")
        .def_static("parse_file", &pygcmc::io::PDBParser::parse_file,
            py::arg("filename"),
            "Parse a PDB file and return the parsed data")
        .def_static("parse_string", &pygcmc::io::PDBParser::parse_string,
            py::arg("pdb_str"),
            "Parse a PDB string and return the parsed data")
        .def_static("parse_to_result", &pygcmc::io::PDBParser::parse_to_result,
            py::arg("filename"),
            py::arg("result"),
            "Parse a PDB file and populate a ParseResult object")
        .def_static("parse_string_to_result", &pygcmc::io::PDBParser::parse_string_to_result,
            py::arg("pdb_str"),
            py::arg("result"),
            "Parse a PDB string and populate a ParseResult object");
    io.attr("PDBParser") = parser;

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
    io.attr("PSFParser") = psf_parser;

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
    io.attr("TOPParser") = top_parser;
}

} // namespace bindings
} // namespace pygcmc


