// src/bindings/model/ModelStructureBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "model/ModelModule.hpp"

namespace py = pybind11;
using namespace pygcmc::model;

namespace pygcmc {
namespace bindings {
namespace model {

void init_structure_bindings(py::module& m, py::module&) {
    // Bind Structure::SecondaryStructure
    auto secondary_structure = py::class_<::pygcmc::model::Structure::SecondaryStructure>(m, "SecondaryStructure")
        .def(py::init<>())
        .def_readwrite("id", &::pygcmc::model::Structure::SecondaryStructure::id)
        .def_readwrite("initResName", &::pygcmc::model::Structure::SecondaryStructure::initResName)
        .def_readwrite("initChainId", &::pygcmc::model::Structure::SecondaryStructure::initChainId)
        .def_readwrite("initSeqNum", &::pygcmc::model::Structure::SecondaryStructure::initSeqNum)
        .def_readwrite("initICode", &::pygcmc::model::Structure::SecondaryStructure::initICode)
        .def_readwrite("endResName", &::pygcmc::model::Structure::SecondaryStructure::endResName)
        .def_readwrite("endChainId", &::pygcmc::model::Structure::SecondaryStructure::endChainId)
        .def_readwrite("endSeqNum", &::pygcmc::model::Structure::SecondaryStructure::endSeqNum)
        .def_readwrite("endICode", &::pygcmc::model::Structure::SecondaryStructure::endICode)
        .def_readwrite("helixClass", &::pygcmc::model::Structure::SecondaryStructure::structureClass);

    // Bind Structure::TerminalInfo
    auto terminal_info = py::class_<::pygcmc::model::Structure::TerminalInfo>(m, "TerminalInfo")
        .def(py::init<>())
        .def_readwrite("chainId", &::pygcmc::model::Structure::TerminalInfo::chainId)
        .def_readwrite("resSeq", &::pygcmc::model::Structure::TerminalInfo::resSeq)
        .def_readwrite("iCode", &::pygcmc::model::Structure::TerminalInfo::iCode)
        .def_readwrite("resName", &::pygcmc::model::Structure::TerminalInfo::resName);

    // Bind Structure to main module
    auto structure = py::class_<::pygcmc::model::Structure>(m, "Structure")
        .def(py::init<>())
        // Directly expose internal members as properties
        .def_property_readonly("atoms", [](const ::pygcmc::model::Structure& s) { return s.get_atoms(); },
            "List of atoms in the structure")
        .def_property_readonly("residues", [](const ::pygcmc::model::Structure& s) { return s.get_residues(); },
            "List of residues in the structure")
        .def_property_readonly("terminals", [](const ::pygcmc::model::Structure& s) { return s.get_terminals(); },
            "List of terminal records")
        .def_property_readonly("helices", [](const ::pygcmc::model::Structure& s) { return s.get_helices(); },
            "Map of chain IDs to helix information")
        .def_property_readonly("sheets", [](const ::pygcmc::model::Structure& s) { return s.get_sheets(); },
            "Map of chain IDs to sheet information")
        .def_property_readonly("ssbonds", [](const ::pygcmc::model::Structure& s) { return s.get_ssbonds(); },
            "List of disulfide bonds")
        .def_property_readonly("box_dimensions", [](const ::pygcmc::model::Structure& s) { return s.get_box_dimensions(); },
            "Box dimensions and angles")
        // Preserve original methods
        .def("add_atom", &::pygcmc::model::Structure::add_atom)
        .def("add_residue", &::pygcmc::model::Structure::add_residue)
        .def("add_terminal", &::pygcmc::model::Structure::add_terminal)
        .def("add_helix", &::pygcmc::model::Structure::add_helix)
        .def("add_sheet", &::pygcmc::model::Structure::add_sheet)
        .def("add_ssbond", &::pygcmc::model::Structure::add_ssbond)
        .def("set_box_dimensions", &::pygcmc::model::Structure::set_box_dimensions)
        .def("clear", &::pygcmc::model::Structure::clear);
}

} // namespace model
} // namespace bindings
} // namespace pygcmc