// src/bindings/model/ModelForcefieldBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "model/ModelModule.hpp"

namespace py = pybind11;
using namespace pygcmc::model;

namespace pygcmc {
namespace bindings {
namespace model {

void init_forcefield_bindings(py::module& m, py::module&) {
    // NonbondedParams
    py::class_<::pygcmc::model::NonbondedParams>(m, "NonbondedParams")
        .def(py::init<>())
        .def_readwrite("nbxmod", &::pygcmc::model::NonbondedParams::nbxmod)
        .def_readwrite("cdiel", &::pygcmc::model::NonbondedParams::cdiel)
        .def_readwrite("fshift", &::pygcmc::model::NonbondedParams::fshift)
        .def_readwrite("vatom", &::pygcmc::model::NonbondedParams::vatom)
        .def_readwrite("vdistance", &::pygcmc::model::NonbondedParams::vdistance)
        .def_readwrite("vfswitch", &::pygcmc::model::NonbondedParams::vfswitch)
        .def_readwrite("cutnb", &::pygcmc::model::NonbondedParams::cutnb)
        .def_readwrite("ctofnb", &::pygcmc::model::NonbondedParams::ctofnb)
        .def_readwrite("ctonnb", &::pygcmc::model::NonbondedParams::ctonnb)
        .def_readwrite("eps", &::pygcmc::model::NonbondedParams::eps)
        .def_readwrite("e14fac", &::pygcmc::model::NonbondedParams::e14fac)
        .def_readwrite("wmin", &::pygcmc::model::NonbondedParams::wmin);

    // LJParams
    py::class_<::pygcmc::model::LJParams>(m, "LJParams")
        .def(py::init<>())
        .def_readwrite("epsilon", &::pygcmc::model::LJParams::epsilon)
        .def_readwrite("rmin_half", &::pygcmc::model::LJParams::rmin_half);

    // BondParams
    py::class_<::pygcmc::model::BondParams>(m, "BondParams")
        .def(py::init<>())
        .def_readwrite("kb", &::pygcmc::model::BondParams::kb)
        .def_readwrite("b0", &::pygcmc::model::BondParams::b0);

    // AngleParams
    py::class_<::pygcmc::model::AngleParams>(m, "AngleParams")
        .def(py::init<>())
        .def_readwrite("ktheta", &::pygcmc::model::AngleParams::ktheta)
        .def_readwrite("theta0", &::pygcmc::model::AngleParams::theta0)
        .def_readwrite("kub", &::pygcmc::model::AngleParams::kub)
        .def_readwrite("s0", &::pygcmc::model::AngleParams::s0);

    // DihedralParams
    py::class_<::pygcmc::model::DihedralParams>(m, "DihedralParams")
        .def(py::init<>())
        .def_readwrite("kchi", &::pygcmc::model::DihedralParams::kchi)
        .def_readwrite("n", &::pygcmc::model::DihedralParams::n)
        .def_readwrite("delta", &::pygcmc::model::DihedralParams::delta);

    // ImproperParams
    py::class_<::pygcmc::model::ImproperParams>(m, "ImproperParams")
        .def(py::init<>())
        .def_readwrite("kpsi", &::pygcmc::model::ImproperParams::kpsi)
        .def_readwrite("psi0", &::pygcmc::model::ImproperParams::psi0);

    // NBFIXParams
    py::class_<::pygcmc::model::NBFIXParams>(m, "NBFIXParams")
        .def(py::init<>())
        .def_readwrite("epsilon", &::pygcmc::model::NBFIXParams::epsilon)
        .def_readwrite("rmin", &::pygcmc::model::NBFIXParams::rmin);

    // Bind ForceField
    py::class_<::pygcmc::model::ForceField>(m, "ForceField")
        .def(py::init<>())
        .def("add_atom_mass", &::pygcmc::model::ForceField::add_atom_mass)
        .def("add_lj_params", &::pygcmc::model::ForceField::add_lj_params)
        .def("add_nbfix", &::pygcmc::model::ForceField::add_nbfix)
        .def("add_bond_params", &::pygcmc::model::ForceField::add_bond_params)
        .def("add_angle_params", &::pygcmc::model::ForceField::add_angle_params)
        .def("add_dihedral_params", &::pygcmc::model::ForceField::add_dihedral_params)
        .def("add_improper_params", &::pygcmc::model::ForceField::add_improper_params)
        .def("get_atom_mass", &::pygcmc::model::ForceField::get_atom_mass)
        .def("get_lj_params", static_cast<const ::pygcmc::model::LJParams& (::pygcmc::model::ForceField::*)(const std::string&) const>(&::pygcmc::model::ForceField::get_lj_params))
        .def("get_nbfix", [](const ::pygcmc::model::ForceField& ff, const std::string& type1, const std::string& type2) {
            auto result = ff.get_nbfix(type1, type2);
            return std::make_tuple(result.first.epsilon, result.first.rmin, result.second);
        }, "Get NBFIX parameters for a pair of atom types. Returns (epsilon, rmin, found)")
        .def("get_bond_params", static_cast<const ::pygcmc::model::BondParams& (::pygcmc::model::ForceField::*)(const std::string&, const std::string&) const>(&::pygcmc::model::ForceField::get_bond_params))
        .def("get_angle_params", static_cast<const ::pygcmc::model::AngleParams& (::pygcmc::model::ForceField::*)(const std::string&, const std::string&, const std::string&) const>(&::pygcmc::model::ForceField::get_angle_params))
        .def("get_dihedral_params", static_cast<const std::vector<::pygcmc::model::DihedralParams>& (::pygcmc::model::ForceField::*)(const std::string&, const std::string&, const std::string&, const std::string&) const>(&::pygcmc::model::ForceField::get_dihedral_params))
        .def("get_improper_params", static_cast<const ::pygcmc::model::ImproperParams& (::pygcmc::model::ForceField::*)(const std::string&, const std::string&, const std::string&, const std::string&) const>(&::pygcmc::model::ForceField::get_improper_params))
        .def("has_atom_mass", &::pygcmc::model::ForceField::has_atom_mass)
        .def("has_lj_params", &::pygcmc::model::ForceField::has_lj_params)
        .def("has_nbfix", &::pygcmc::model::ForceField::has_nbfix)
        .def("has_bond_params", &::pygcmc::model::ForceField::has_bond_params)
        .def("has_angle_params", &::pygcmc::model::ForceField::has_angle_params)
        .def("has_dihedral_params", &::pygcmc::model::ForceField::has_dihedral_params)
        .def("has_improper_params", &::pygcmc::model::ForceField::has_improper_params)
        .def("get_num_atom_types", &::pygcmc::model::ForceField::get_num_atom_types)
        .def("get_num_lj_params", &::pygcmc::model::ForceField::get_num_lj_params)
        .def("get_num_nbfix", &::pygcmc::model::ForceField::get_num_nbfix)
        .def("get_num_bond_types", &::pygcmc::model::ForceField::get_num_bond_types)
        .def("get_num_angle_types", &::pygcmc::model::ForceField::get_num_angle_types)
        .def("get_num_dihedral_types", &::pygcmc::model::ForceField::get_num_dihedral_types)
        .def("get_num_improper_types", &::pygcmc::model::ForceField::get_num_improper_types)
        .def("get_nonbonded_params", static_cast<const ::pygcmc::model::NonbondedParams& (::pygcmc::model::ForceField::*)() const>(&::pygcmc::model::ForceField::get_nonbonded_params))
        .def_static("makeTypePair", &::pygcmc::model::ForceField::makeTypePair)
        .def_static("makeTypeTriple", &::pygcmc::model::ForceField::makeTypeTriple)
        .def_static("makeTypeQuad", &::pygcmc::model::ForceField::makeTypeQuad)
        // Property accessors - using new _map methods for accessing full parameter maps
        .def_property_readonly("atom_masses", static_cast<const std::map<std::string, double>& (::pygcmc::model::ForceField::*)() const>(&::pygcmc::model::ForceField::get_atom_masses))
        .def_property_readonly("lj_params", static_cast<const std::map<std::string, ::pygcmc::model::LJParams>& (::pygcmc::model::ForceField::*)() const>(&::pygcmc::model::ForceField::get_lj_params_map))
        .def_property_readonly("nbfix", static_cast<const std::map<std::pair<std::string, std::string>, ::pygcmc::model::NBFIXParams>& (::pygcmc::model::ForceField::*)() const>(&::pygcmc::model::ForceField::get_nbfix_map))
        .def_property_readonly("bond_params", static_cast<const std::map<std::pair<std::string, std::string>, ::pygcmc::model::BondParams>& (::pygcmc::model::ForceField::*)() const>(&::pygcmc::model::ForceField::get_bond_params_map))
        .def_property_readonly("angle_params", static_cast<const std::map<std::tuple<std::string, std::string, std::string>, ::pygcmc::model::AngleParams>& (::pygcmc::model::ForceField::*)() const>(&::pygcmc::model::ForceField::get_angle_params_map))
        .def_property_readonly("dihedral_params", static_cast<const std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<::pygcmc::model::DihedralParams>>& (::pygcmc::model::ForceField::*)() const>(&::pygcmc::model::ForceField::get_dihedral_params_map))
        .def_property_readonly("improper_params", static_cast<const std::map<std::tuple<std::string, std::string, std::string, std::string>, ::pygcmc::model::ImproperParams>& (::pygcmc::model::ForceField::*)() const>(&::pygcmc::model::ForceField::get_improper_params_map));
}

} // namespace model
} // namespace bindings
} // namespace pygcmc