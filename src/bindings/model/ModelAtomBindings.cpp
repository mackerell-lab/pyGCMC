// src/bindings/model/ModelAtomBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "model/ModelModule.hpp"

namespace py = pybind11;
using namespace pygcmc::model;
using namespace pygcmc::model::atom;

namespace pygcmc {
namespace bindings {
namespace model {

void init_atom_bindings(py::module& m, py::module& model_module) {
    // ----------------------------------------------------------
    // Internal base class bindings (no Python exposure)
    // Needed so that derived classes can list them as bases
    // ----------------------------------------------------------
    py::class_<::pygcmc::model::atom::AtomCore, std::shared_ptr<::pygcmc::model::atom::AtomCore>>(m, "_AtomCore");

    // Bind Atom class with base AtomCore to enable inherited method bindings
    py::class_<::pygcmc::model::Atom, ::pygcmc::model::atom::AtomCore, std::shared_ptr<::pygcmc::model::Atom>>(model_module, "Atom")
        .def(py::init<>())
        .def("get_bynu", &::pygcmc::model::Atom::get_bynu)
        .def("get_type", &::pygcmc::model::Atom::get_type)
        .def("get_resname", &::pygcmc::model::Atom::get_resname)
        .def("get_ires", &::pygcmc::model::Atom::get_ires)
        .def("get_chain", &::pygcmc::model::Atom::get_chain)
        .def("get_coor", &::pygcmc::model::Atom::get_coor)
        .def("get_x", &::pygcmc::model::Atom::get_x)
        .def("get_y", &::pygcmc::model::Atom::get_y)
        .def("get_z", &::pygcmc::model::Atom::get_z)
        .def("is_hetatm", &::pygcmc::model::Atom::is_hetatm)
        .def("get_occupancy", &::pygcmc::model::Atom::get_occupancy)
        .def("get_tempfactor", &::pygcmc::model::Atom::get_tempfactor)
        .def("get_segid", &::pygcmc::model::Atom::get_segid)
        .def("get_inscode", &::pygcmc::model::Atom::get_inscode)
        .def("get_formatted_atom_name", &::pygcmc::model::Atom::get_formatted_atom_name)
        .def("get_name", &::pygcmc::model::Atom::get_name)
        .def("get_residue_id", &::pygcmc::model::Atom::get_residue_id)
        .def("atom_name", &::pygcmc::model::Atom::atom_name)
        .def("residue_name", &::pygcmc::model::Atom::residue_name)
        .def("residue_number", &::pygcmc::model::Atom::residue_number)
        .def("get_mass", &::pygcmc::model::Atom::get_mass)
        .def("get_charge", &::pygcmc::model::Atom::get_charge)
        .def("get_alpha", &::pygcmc::model::Atom::get_alpha)
        .def("get_thole", &::pygcmc::model::Atom::get_thole)
        .def("set_residue_id", &::pygcmc::model::Atom::set_residue_id)
        .def("set_coor", &::pygcmc::model::Atom::set_coor)
        .def("set_mass_charge", &::pygcmc::model::Atom::set_mass_charge)
        .def("set_lj_params", &::pygcmc::model::Atom::set_lj_params)
        .def("set_drude_params", &::pygcmc::model::Atom::set_drude_params)
        .def("set_alpha", &::pygcmc::model::Atom::set_alpha)
        .def("set_thole", &::pygcmc::model::Atom::set_thole)
        .def("set_occupancy", &::pygcmc::model::Atom::set_occupancy)
        .def("set_tempfactor", &::pygcmc::model::Atom::set_tempfactor)
        .def("set_element", &::pygcmc::model::Atom::set_element)
        .def("set_charge_string", &::pygcmc::model::Atom::set_charge_string)
        .def("set_chain", &::pygcmc::model::Atom::set_chain)
        .def("set_hetatm", &::pygcmc::model::Atom::set_hetatm)
        .def("set_bynu", &::pygcmc::model::Atom::set_bynu)
        .def("set_type", &::pygcmc::model::Atom::set_type)
        .def("set_resname", &::pygcmc::model::Atom::set_resname)
        .def("set_ires", &::pygcmc::model::Atom::set_ires)
        .def("set_segid", &::pygcmc::model::Atom::set_segid)
        .def("set_altloc", &::pygcmc::model::Atom::set_altloc)
        .def("set_inscode", &::pygcmc::model::Atom::set_inscode)
        .def("has_lj_params", &::pygcmc::model::Atom::has_lj_params)
        .def("is_valid", &::pygcmc::model::Atom::is_valid);
}

} // namespace model
} // namespace bindings
} // namespace pygcmc
