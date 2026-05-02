// src/bindings/model/ModelResidueBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "model/ModelModule.hpp"

namespace py = pybind11;
using namespace pygcmc::model;

namespace pygcmc {
namespace bindings {
namespace model {

void init_residue_bindings(py::module&, py::module& model_module) {
    // Bind Residue class
    py::class_<::pygcmc::model::Residue, std::shared_ptr<::pygcmc::model::Residue>>(model_module, "Residue")
        .def(py::init<>())
        .def("get_resname", &::pygcmc::model::Residue::get_resname)
        .def("set_resname", &::pygcmc::model::Residue::set_resname)
        .def("get_ires", &::pygcmc::model::Residue::get_ires)
        .def("get_segid", &::pygcmc::model::Residue::get_segid)
        .def("get_iseg", &::pygcmc::model::Residue::get_iseg)
        .def("get_chain", &::pygcmc::model::Residue::get_chain)
        .def("get_inscode", &::pygcmc::model::Residue::get_inscode)
        .def("get_atoms", &::pygcmc::model::Residue::get_atoms)
        .def("add_atom", py::overload_cast<const ::pygcmc::model::Atom&>(&::pygcmc::model::Residue::add_atom))
        .def("add_atom", py::overload_cast<std::shared_ptr<::pygcmc::model::Atom>>(&::pygcmc::model::Residue::add_atom))
        .def("find_atom", &::pygcmc::model::Residue::find_atom)
        .def("atom_count", &::pygcmc::model::Residue::atom_count)
        .def("is_valid", &::pygcmc::model::Residue::is_valid)
        .def("calculate_center_of_mass", &::pygcmc::model::Residue::calculate_center_of_mass)
        .def("get_center_of_mass", &::pygcmc::model::Residue::get_center_of_mass)
        .def("has_atom_type", &::pygcmc::model::Residue::has_atom_type)
        .def("get_residue_id", &::pygcmc::model::Residue::get_residue_id)
        .def("set_residue_id", &::pygcmc::model::Residue::set_residue_id)
        .def("find_atom_by_pdb_name", &::pygcmc::model::Residue::find_atom_by_pdb_name)
        .def("get_atom_range", &::pygcmc::model::Residue::get_atom_range);
}

} // namespace model
} // namespace bindings
} // namespace pygcmc
