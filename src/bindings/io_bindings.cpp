// src/bindings/io_bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "io/pdbParser.hpp"
#include "model/atom.hpp"
#include "model/residue.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {

void init_io(py::module& m) {
    // Create io submodule
    auto io = m.def_submodule("io", "Input/Output operations");
    
    // Bind PDBParser and its ParseResult in both main module and io submodule
    py::class_<io::PDBParser::ParseResult> parse_result(m, "PDBParseResult");
    parse_result
        .def(py::init<>())
        .def_readwrite("atoms", &io::PDBParser::ParseResult::atoms)
        .def_readwrite("residues", &io::PDBParser::ParseResult::residues)
        .def_readwrite("helices", &io::PDBParser::ParseResult::helices)
        .def_readwrite("sheets", &io::PDBParser::ParseResult::sheets)
        .def_readwrite("ssbonds", &io::PDBParser::ParseResult::ssbonds)
        .def_readwrite("boxDimensions", &io::PDBParser::ParseResult::boxDimensions);

    py::class_<io::PDBParser> parser(m, "PDBParser");
    parser
        .def_static("parse_file", &io::PDBParser::parseFile,
            py::arg("filename"),
            "Parse a PDB file and return the parsed data")
        .def_static("parse_string", &io::PDBParser::parseString,
            py::arg("pdb_str"),
            "Parse a PDB string and return the parsed data");

    // Also bind to io submodule
    io.attr("PDBParser") = parser;
    io.attr("PDBParseResult") = parse_result;
}

} // namespace bindings
} // namespace pygcmc


