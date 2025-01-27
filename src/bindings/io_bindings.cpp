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
        .def_static("parse_file", &pygcmc::io::PDBParser::parseFile,
            py::arg("filename"),
            "Parse a PDB file and return the parsed data")
        .def_static("parse_string", &pygcmc::io::PDBParser::parseString,
            py::arg("pdb_str"),
            "Parse a PDB string and return the parsed data");
    io.attr("PDBParser") = parser;

    // PSFParser bindings
    auto psf_parser = py::class_<pygcmc::io::PSFParser>(m, "PSFParser")
        .def(py::init<>())
        .def("parse_to_topology", &pygcmc::io::PSFParser::parse_to_topology);
    io.attr("PSFParser") = psf_parser;

    // TopParser bindings
    auto top_parser = py::class_<pygcmc::io::TopParser>(m, "TopParser")
        .def(py::init<>())
        .def("parse_to_topology", &pygcmc::io::TopParser::parse_to_topology);
    io.attr("TopParser") = top_parser;
}

} // namespace bindings
} // namespace pygcmc


