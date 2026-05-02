// src/bindings/model/ModelMoleculeBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "model/ModelModule.hpp"

namespace py = pybind11;
using namespace pygcmc::model;

namespace pygcmc {
namespace bindings {
namespace model {

void init_molecule_bindings(py::module&, py::module& model_module) {
    // Bind Molecular class
    py::class_<::pygcmc::model::Molecular, std::shared_ptr<::pygcmc::model::Molecular>>(model_module, "Molecular")
        .def(py::init<>())
        .def_readwrite("atoms", &::pygcmc::model::Molecular::atoms)
        .def_readwrite("residues", &::pygcmc::model::Molecular::residues)
        .def_readwrite("terminals", &::pygcmc::model::Molecular::terminals)
        .def_readwrite("helices", &::pygcmc::model::Molecular::helices)
        .def_readwrite("sheets", &::pygcmc::model::Molecular::sheets)
        .def_readwrite("ssbonds", &::pygcmc::model::Molecular::ssbonds)
        .def_readwrite("boxDimensions", &::pygcmc::model::Molecular::boxDimensions)
        .def_readwrite("topology_atoms", &::pygcmc::model::Molecular::topology_atoms)
        .def_readwrite("topology_residues", &::pygcmc::model::Molecular::topology_residues)
        .def_readwrite("segments", &::pygcmc::model::Molecular::segments)
        .def_readwrite("bonds", &::pygcmc::model::Molecular::bonds)
        .def_readwrite("angles", &::pygcmc::model::Molecular::angles)
        .def_readwrite("dihedrals", &::pygcmc::model::Molecular::dihedrals)
        .def_readwrite("donors", &::pygcmc::model::Molecular::donors)
        .def_readwrite("acceptors", &::pygcmc::model::Molecular::acceptors)
        .def_readwrite("exclusions", &::pygcmc::model::Molecular::exclusions)
        .def_readwrite("groups", &::pygcmc::model::Molecular::groups)
        .def_readwrite("cmaps", &::pygcmc::model::Molecular::cmaps)
        .def_readwrite("standard_cmaps", &::pygcmc::model::Molecular::standard_cmaps)
        .def_readwrite("titles", &::pygcmc::model::Molecular::titles)
        .def_readwrite("segment_map", &::pygcmc::model::Molecular::segment_map)
        .def_readwrite("residue_map", &::pygcmc::model::Molecular::residue_map)
        .def_readwrite("atom_map", &::pygcmc::model::Molecular::atom_map)
        .def("get_num_atoms", &::pygcmc::model::Molecular::get_num_atoms)
        .def("get_num_residues", &::pygcmc::model::Molecular::get_num_residues)
        .def("get_num_segments", &::pygcmc::model::Molecular::get_num_segments)
        .def("get_num_bonds", &::pygcmc::model::Molecular::get_num_bonds)
        .def("get_num_angles", &::pygcmc::model::Molecular::get_num_angles)
        .def("get_num_dihedrals", &::pygcmc::model::Molecular::get_num_dihedrals)
        .def("get_num_impropers", &::pygcmc::model::Molecular::get_num_impropers)
        .def("get_num_standard_cmaps", &::pygcmc::model::Molecular::get_num_standard_cmaps)
        .def("add_standard_cmap", &::pygcmc::model::Molecular::add_standard_cmap)
        .def("clear", &::pygcmc::model::Molecular::clear);
}

} // namespace model
} // namespace bindings
} // namespace pygcmc
