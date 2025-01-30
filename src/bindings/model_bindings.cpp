// src/bindings/model_bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "model/atom.hpp"
#include "model/residue.hpp"
#include "model/topology.hpp"
#include "../model/forcefield.hpp"

namespace py = pybind11;
using namespace pygcmc;

namespace pygcmc {
namespace bindings {

void init_model(py::module& m) {
    // Create model submodule
    auto model = m.def_submodule("model", "Data model classes");
    
    // Bind Atom class
    py::class_<model::Atom, std::shared_ptr<model::Atom>>(model, "PDBAtom")
        .def(py::init<>())
        .def("getBynu", &model::Atom::getBynu)
        .def("getType", &model::Atom::getType)
        .def("getResname", &model::Atom::getResname)
        .def("getIres", &model::Atom::getIres)
        .def("getChain", &model::Atom::getChain)
        .def("getCoor", &model::Atom::getCoor)
        .def("isHetatm", &model::Atom::isHetatm)
        .def("getOccupancy", &model::Atom::getOccupancy)
        .def("getTempfactor", &model::Atom::getTempfactor)
        .def("getSegid", &model::Atom::getSegid)
        .def("getInscode", &model::Atom::getInscode)
        .def("getFormattedAtomName", &model::Atom::getFormattedAtomName)
        .def("getResidueID", &model::Atom::getResidueID)
        .def("getMass", &model::Atom::getMass)
        .def("setResidueID", &model::Atom::setResidueID)
        .def("setCoor", &model::Atom::setCoor)
        .def("setMassCharge", &model::Atom::setMassCharge)
        .def("setLJParams", &model::Atom::setLJParams)
        .def("setOccupancy", &model::Atom::setOccupancy)
        .def("setTempfactor", &model::Atom::setTempfactor)
        .def("setElement", &model::Atom::setElement)
        .def("setChargeString", &model::Atom::setChargeString)
        .def("setChain", &model::Atom::setChain)
        .def("setHetatm", &model::Atom::setHetatm)
        .def("setBynu", &model::Atom::setBynu)
        .def("setType", &model::Atom::setType)
        .def("setResname", &model::Atom::setResname)
        .def("setIres", &model::Atom::setIres)
        .def("setSegid", &model::Atom::setSegid)
        .def("setAltloc", &model::Atom::setAltloc)
        .def("setInscode", &model::Atom::setInscode)
        .def("hasLJParams", &model::Atom::hasLJParams)
        .def("isValid", &model::Atom::isValid);

    // Bind PDB Residue class
    py::class_<model::Residue, std::shared_ptr<model::Residue>>(model, "PDBResidue")
        .def(py::init<>())
        .def("getResname", &model::Residue::getResname)
        .def("getIres", &model::Residue::getIres)
        .def("getSegid", &model::Residue::getSegid)
        .def("getIseg", &model::Residue::getIseg)
        .def("getChain", &model::Residue::getChain)
        .def("getInscode", &model::Residue::getInscode)
        .def("getAtoms", &model::Residue::getAtoms)
        .def("addAtom", py::overload_cast<const model::Atom&>(&model::Residue::addAtom))
        .def("addAtom", py::overload_cast<std::shared_ptr<model::Atom>>(&model::Residue::addAtom))
        .def("findAtom", &model::Residue::findAtom)
        .def("atomCount", &model::Residue::atomCount)
        .def("isValid", &model::Residue::isValid)
        .def("calculateCenterOfMass", &model::Residue::calculateCenterOfMass)
        .def("getCenterOfMass", &model::Residue::getCenterOfMass)
        .def("hasAtomType", &model::Residue::hasAtomType)
        .def("getResidueID", &model::Residue::getResidueID)
        .def("setResidueID", &model::Residue::setResidueID)
        .def("findAtomByPDBName", &model::Residue::findAtomByPDBName)
        .def("updateAtomMap", &model::Residue::updateAtomMap)
        .def("getAtomRange", &model::Residue::getAtomRange);

    // Bind Topology class
    py::class_<model::Topology>(model, "Topology")
        .def(py::init<>())
        .def("add_atom", &model::Topology::add_atom)
        .def("add_bond", &model::Topology::add_bond)
        .def("add_angle", &model::Topology::add_angle)
        .def("add_dihedral", &model::Topology::add_dihedral)
        .def("add_improper", &model::Topology::add_improper)
        .def("add_donor", &model::Topology::add_donor)
        .def("add_acceptor", &model::Topology::add_acceptor)
        .def("add_nonbonded_exclusion", &model::Topology::add_nonbonded_exclusion)
        .def("add_group", &model::Topology::add_group)
        .def("add_cmap", static_cast<void (model::Topology::*)(const std::array<int, 8>&)>(&model::Topology::add_cmap), "Add CMAP in CHARMM format (8 atoms)")
        .def("add_cmap", static_cast<void (model::Topology::*)(const std::array<int, 5>&, int)>(&model::Topology::add_cmap), "Add CMAP in GROMACS format (5 atoms)", py::arg("atoms"), py::arg("function_type") = 1)
        .def("add_title", &model::Topology::add_title)
        .def("get_num_atoms", &model::Topology::get_num_atoms)
        .def("get_num_residues", &model::Topology::get_num_residues)
        .def("get_num_segments", &model::Topology::get_num_segments)
        .def("get_num_bonds", &model::Topology::get_num_bonds)
        .def("get_num_angles", &model::Topology::get_num_angles)
        .def("get_num_dihedrals", &model::Topology::get_num_dihedrals)
        .def("get_num_impropers", &model::Topology::get_num_impropers)
        .def("get_num_donors", &model::Topology::get_num_donors)
        .def("get_num_acceptors", &model::Topology::get_num_acceptors)
        .def("get_num_cmaps", &model::Topology::get_num_cmaps)
        .def("get_num_groups", &model::Topology::get_num_groups)
        .def("get_group", &model::Topology::get_group, py::return_value_policy::reference_internal)
        .def("has_atom", &model::Topology::has_atom)
        .def("has_group", &model::Topology::has_group)
        .def("has_bond", &model::Topology::has_bond)
        .def("has_angle", &model::Topology::has_angle)
        .def("has_dihedral", &model::Topology::has_dihedral)
        .def("has_improper", &model::Topology::has_improper)
        .def("has_donor", &model::Topology::has_donor)
        .def("has_acceptor", &model::Topology::has_acceptor)
        .def("has_cmap", &model::Topology::has_cmap)
        .def("get_residue", &model::Topology::get_residue, py::return_value_policy::reference_internal)
        .def("get_atom", &model::Topology::get_atom, py::return_value_policy::reference_internal)
        .def("find_residue", &model::Topology::find_residue);

    // Bind TopologyResidue
    py::class_<model::TopologyResidue>(model, "TopologyResidue")
        .def_readonly("name", &model::TopologyResidue::name)
        .def_readonly("number", &model::TopologyResidue::number)
        .def_readonly("atoms", &model::TopologyResidue::atoms)
        .def_readonly("segment", &model::TopologyResidue::segment);

    // Bind TopologyAtom
    py::class_<model::TopologyAtom>(model, "TopologyAtom")
        .def_readonly("name", &model::TopologyAtom::name)
        .def_readonly("type", &model::TopologyAtom::type)
        .def_readonly("charge", &model::TopologyAtom::charge)
        .def_readonly("mass", &model::TopologyAtom::mass);

    // Bind TopologyGroup
    py::class_<model::TopologyGroup>(model, "TopologyGroup")
        .def_readonly("id", &model::TopologyGroup::id)
        .def_readonly("atoms", &model::TopologyGroup::atoms)
        .def_readonly("type", &model::TopologyGroup::type);

    // NonbondedParams
    py::class_<NonbondedParams>(m, "NonbondedParams")
        .def(py::init<>())
        .def_readwrite("nbxmod", &NonbondedParams::nbxmod)
        .def_readwrite("cdiel", &NonbondedParams::cdiel)
        .def_readwrite("fshift", &NonbondedParams::fshift)
        .def_readwrite("vatom", &NonbondedParams::vatom)
        .def_readwrite("vdistance", &NonbondedParams::vdistance)
        .def_readwrite("vfswitch", &NonbondedParams::vfswitch)
        .def_readwrite("cutnb", &NonbondedParams::cutnb)
        .def_readwrite("ctofnb", &NonbondedParams::ctofnb)
        .def_readwrite("ctonnb", &NonbondedParams::ctonnb)
        .def_readwrite("eps", &NonbondedParams::eps)
        .def_readwrite("e14fac", &NonbondedParams::e14fac)
        .def_readwrite("wmin", &NonbondedParams::wmin);

    // LJParams
    py::class_<LJParams>(m, "LJParams")
        .def(py::init<>())
        .def_readwrite("epsilon", &LJParams::epsilon)
        .def_readwrite("rmin", &LJParams::rmin);

    // BondParams
    py::class_<BondParams>(m, "BondParams")
        .def(py::init<>())
        .def_readwrite("kb", &BondParams::kb)
        .def_readwrite("b0", &BondParams::b0);

    // AngleParams
    py::class_<AngleParams>(m, "AngleParams")
        .def(py::init<>())
        .def_readwrite("ktheta", &AngleParams::ktheta)
        .def_readwrite("theta0", &AngleParams::theta0)
        .def_readwrite("kub", &AngleParams::kub)
        .def_readwrite("s0", &AngleParams::s0);

    // DihedralParams
    py::class_<DihedralParams>(m, "DihedralParams")
        .def(py::init<>())
        .def_readwrite("kchi", &DihedralParams::kchi)
        .def_readwrite("n", &DihedralParams::n)
        .def_readwrite("delta", &DihedralParams::delta);

    // ImproperParams
    py::class_<ImproperParams>(m, "ImproperParams")
        .def(py::init<>())
        .def_readwrite("kpsi", &ImproperParams::kpsi)
        .def_readwrite("psi0", &ImproperParams::psi0);

    // ForceField
    py::class_<ForceField>(m, "ForceField")
        .def(py::init<>())
        .def("add_atom_mass", &ForceField::add_atom_mass)
        .def("add_lj_params", &ForceField::add_lj_params)
        .def("add_nbfix", &ForceField::add_nbfix)
        .def("add_bond_params", &ForceField::add_bond_params)
        .def("add_angle_params", &ForceField::add_angle_params)
        .def("add_dihedral_params", &ForceField::add_dihedral_params)
        .def("add_improper_params", &ForceField::add_improper_params)
        .def("get_atom_mass", &ForceField::get_atom_mass)
        .def("get_lj_params", static_cast<const LJParams& (ForceField::*)(const std::string&) const>(&ForceField::get_lj_params))
        .def("get_nbfix", static_cast<std::pair<double, bool> (ForceField::*)(const std::string&, const std::string&) const>(&ForceField::get_nbfix))
        .def("get_bond_params", static_cast<const BondParams& (ForceField::*)(const std::string&, const std::string&) const>(&ForceField::get_bond_params))
        .def("get_angle_params", static_cast<const AngleParams& (ForceField::*)(const std::string&, const std::string&, const std::string&) const>(&ForceField::get_angle_params))
        .def("get_dihedral_params", static_cast<const std::vector<DihedralParams>& (ForceField::*)(const std::string&, const std::string&, const std::string&, const std::string&) const>(&ForceField::get_dihedral_params))
        .def("get_improper_params", static_cast<const ImproperParams& (ForceField::*)(const std::string&, const std::string&, const std::string&, const std::string&) const>(&ForceField::get_improper_params))
        .def("has_atom_mass", &ForceField::has_atom_mass)
        .def("has_lj_params", &ForceField::has_lj_params)
        .def("has_nbfix", &ForceField::has_nbfix)
        .def("has_bond_params", &ForceField::has_bond_params)
        .def("has_angle_params", &ForceField::has_angle_params)
        .def("has_dihedral_params", &ForceField::has_dihedral_params)
        .def("has_improper_params", &ForceField::has_improper_params)
        .def("get_num_atom_types", &ForceField::get_num_atom_types)
        .def("get_num_lj_params", &ForceField::get_num_lj_params)
        .def("get_num_nbfix", &ForceField::get_num_nbfix)
        .def("get_num_bond_types", &ForceField::get_num_bond_types)
        .def("get_num_angle_types", &ForceField::get_num_angle_types)
        .def("get_num_dihedral_types", &ForceField::get_num_dihedral_types)
        .def("get_num_improper_types", &ForceField::get_num_improper_types)
        .def("get_nonbonded_params", static_cast<const NonbondedParams& (ForceField::*)() const>(&ForceField::get_nonbonded_params))
        .def_static("makeTypePair", &ForceField::makeTypePair)
        .def_static("makeTypeTriple", &ForceField::makeTypeTriple)
        .def_static("makeTypeQuad", &ForceField::makeTypeQuad)
        // Add property accessors for the maps
        .def_property_readonly("atom_masses", static_cast<const std::map<std::string, double>& (ForceField::*)() const>(&ForceField::get_atom_masses))
        .def_property_readonly("lj_params", static_cast<const std::map<std::string, LJParams>& (ForceField::*)() const>(&ForceField::get_lj_params))
        .def_property_readonly("nbfix", static_cast<const std::map<std::pair<std::string, std::string>, double>& (ForceField::*)() const>(&ForceField::get_nbfix))
        .def_property_readonly("bond_params", static_cast<const std::map<std::pair<std::string, std::string>, BondParams>& (ForceField::*)() const>(&ForceField::get_bond_params))
        .def_property_readonly("angle_params", static_cast<const std::map<std::tuple<std::string, std::string, std::string>, AngleParams>& (ForceField::*)() const>(&ForceField::get_angle_params))
        .def_property_readonly("dihedral_params", static_cast<const std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<DihedralParams>>& (ForceField::*)() const>(&ForceField::get_dihedral_params))
        .def_property_readonly("improper_params", static_cast<const std::map<std::tuple<std::string, std::string, std::string, std::string>, ImproperParams>& (ForceField::*)() const>(&ForceField::get_improper_params));
}

} // namespace bindings
} // namespace pygcmc 