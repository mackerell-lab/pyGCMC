// modules/bindings/src/parser_bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "pygcmc/core/io/pdb_parser.hpp"
#include "pygcmc/core/io/top_parser.hpp"
#include "pygcmc/core/io/ff_parser.hpp"
#include "pygcmc/core/project_atom.hpp"

namespace py = pybind11;
using namespace pygcmc::core;

void init_parser_bindings(py::module& m) {
    // Bind PDBParser class
    py::class_<io::PDBParser>(m, "PDBParser")
        .def(py::init<>())
        .def_static("parse", &io::PDBParser::parse);

    // Bind Residue struct for PDB with proper atom handling
    py::class_<io::Residue>(m, "PDBResidue")
        .def(py::init<>())
        .def_readwrite("name", &io::Residue::name)
        .def_readwrite("sequence_number", &io::Residue::sequence_number)
        .def_readwrite("chain_id", &io::Residue::chain_id)
        .def_readwrite("atoms", &io::Residue::atoms);

    // Bind TopParser class
    py::class_<io::TopParser>(m, "TopParser")
        .def(py::init<>())
        .def("parse", &io::TopParser::parse)
        .def("parse_with_includes", &io::TopParser::parse_with_includes)
        .def("update_pdb_atoms", [](io::TopParser& self, py::list atoms) -> int {
            std::vector<io::PDBAtom*> c_atoms;
            c_atoms.reserve(py::len(atoms));
            
            for(auto item : atoms) {
                try {
                    // Try as raw PDBAtom first
                    auto& atom = item.cast<io::PDBAtom&>();
                    c_atoms.push_back(&atom);
                } catch (const py::cast_error&) {
                    try {
                        // Try as ProjectAtom
                        auto& wrapper = item.cast<ProjectAtom&>();
                        c_atoms.push_back(wrapper.get_ptr().get());
                    } catch (const py::cast_error&) {
                        // Try as shared_ptr
                        auto shared_atom = item.cast<std::shared_ptr<io::PDBAtom>>();
                        c_atoms.push_back(shared_atom.get());
                    }
                }
            }
            return self.update_pdb_atoms(c_atoms);
        }, py::arg("atoms"));

    // Bind FFParser class
    py::class_<io::FFParser>(m, "FFParser")
        .def(py::init<>())
        .def("parse", &io::FFParser::parse)
        .def("get_nonbonded_params", &io::FFParser::get_nonbonded_params)
        .def("get_nbfix_params", &io::FFParser::get_nbfix_params)
        .def("update_pdb_atoms", [](io::FFParser& self, py::list atoms) -> int {
            std::vector<io::PDBAtom*> atom_ptrs;
            atom_ptrs.reserve(py::len(atoms));
            
            for (auto item : atoms) {
                try {
                    // Try as raw PDBAtom first
                    auto& atom = item.cast<io::PDBAtom&>();
                    atom_ptrs.push_back(&atom);
                } catch (const py::cast_error&) {
                    try {
                        // Try as ProjectAtom
                        auto& wrapper = item.cast<ProjectAtom&>();
                        atom_ptrs.push_back(wrapper.get_ptr().get());
                    } catch (const py::cast_error&) {
                        // Try as shared_ptr
                        auto shared_atom = item.cast<std::shared_ptr<io::PDBAtom>>();
                        atom_ptrs.push_back(shared_atom.get());
                    }
                }
            }
            return self.update_pdb_atoms(atom_ptrs);
        })
        // Add global nonbonded parameter getters
        .def("get_cutnb", &io::FFParser::get_cutnb)
        .def("get_ctofnb", &io::FFParser::get_ctofnb)
        .def("get_ctonnb", &io::FFParser::get_ctonnb)
        .def("get_eps", &io::FFParser::get_eps)
        .def("get_e14fac", &io::FFParser::get_e14fac)
        .def("get_wmin", &io::FFParser::get_wmin)
        // Add static method for merging NBFIX parameters
        .def_static("merge_nbfix_params", [](py::list parsers) {
            std::vector<const io::FFParser*> parser_ptrs;
            for (auto item : parsers) {
                parser_ptrs.push_back(item.cast<io::FFParser*>());
            }
            return io::FFParser::merge_nbfix_params(parser_ptrs);
        });
} 