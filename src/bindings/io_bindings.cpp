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
    
    // Bind PDBParser
    py::class_<io::PDBParser::ParseResult>(io, "PDBParseResult")
        .def(py::init<>())
        .def_readwrite("atoms", &io::PDBParser::ParseResult::atoms)
        .def_readwrite("residues", &io::PDBParser::ParseResult::residues)
        .def_readwrite("helices", &io::PDBParser::ParseResult::helices)
        .def_readwrite("sheets", &io::PDBParser::ParseResult::sheets)
        .def_readwrite("ssbonds", &io::PDBParser::ParseResult::ssbonds);

    py::class_<io::PDBParser>(io, "PDBParser")
        .def_static("parse_file", &io::PDBParser::parseFile,
            py::arg("filename"),
            "Parse a PDB file and return the parsed data")
        .def_static("parse_string", &io::PDBParser::parseString,
            py::arg("pdb_str"),
            "Parse a PDB string and return the parsed data");
}

} // namespace bindings
} // namespace pygcmc


