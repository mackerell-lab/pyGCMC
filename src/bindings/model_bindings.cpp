// src/bindings/model_bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "model/atom.hpp"
#include "model/residue.hpp"
#include "model/topology.hpp"
#include "model/structure.hpp"
#include "model/forcefield.hpp"

namespace py = pybind11;
// using namespace pygcmc;

namespace pygcmc {
namespace bindings {

void init_model(py::module& m) {
    // Create model submodule
    auto model = m.def_submodule("model", "Data model classes");
    
    // Bind Structure::SecondaryStructure
    auto secondary_structure = py::class_<model::Structure::SecondaryStructure>(m, "SecondaryStructure")
        .def(py::init<>())
        .def_readwrite("id", &model::Structure::SecondaryStructure::id)
        .def_readwrite("initResName", &model::Structure::SecondaryStructure::initResName)
        .def_readwrite("initChainId", &model::Structure::SecondaryStructure::initChainId)
        .def_readwrite("initSeqNum", &model::Structure::SecondaryStructure::initSeqNum)
        .def_readwrite("initICode", &model::Structure::SecondaryStructure::initICode)
        .def_readwrite("endResName", &model::Structure::SecondaryStructure::endResName)
        .def_readwrite("endChainId", &model::Structure::SecondaryStructure::endChainId)
        .def_readwrite("endSeqNum", &model::Structure::SecondaryStructure::endSeqNum)
        .def_readwrite("endICode", &model::Structure::SecondaryStructure::endICode)
        .def_readwrite("helixClass", &model::Structure::SecondaryStructure::structureClass);

    // Bind Structure::TerminalInfo
    auto terminal_info = py::class_<model::Structure::TerminalInfo>(m, "TerminalInfo")
        .def(py::init<>())
        .def_readwrite("chainId", &model::Structure::TerminalInfo::chainId)
        .def_readwrite("resSeq", &model::Structure::TerminalInfo::resSeq)
        .def_readwrite("iCode", &model::Structure::TerminalInfo::iCode)
        .def_readwrite("resName", &model::Structure::TerminalInfo::resName);

    // Bind Structure to main module
    auto structure = py::class_<model::Structure>(m, "Structure")
        .def(py::init<>())
        // 直接暴露内部成员作为属性
        .def_property_readonly("atoms", [](const model::Structure& s) { return s.get_atoms(); },
            "List of atoms in the structure")
        .def_property_readonly("residues", [](const model::Structure& s) { return s.get_residues(); },
            "List of residues in the structure")
        .def_property_readonly("terminals", [](const model::Structure& s) { return s.get_terminals(); },
            "List of terminal records")
        .def_property_readonly("helices", [](const model::Structure& s) { return s.get_helices(); },
            "Map of chain IDs to helix information")
        .def_property_readonly("sheets", [](const model::Structure& s) { return s.get_sheets(); },
            "Map of chain IDs to sheet information")
        .def_property_readonly("ssbonds", [](const model::Structure& s) { return s.get_ssbonds(); },
            "List of disulfide bonds")
        .def_property_readonly("box_dimensions", [](const model::Structure& s) { return s.get_box_dimensions(); },
            "Box dimensions and angles")
        // 保留原有方法
        .def("add_atom", &model::Structure::add_atom)
        .def("add_residue", &model::Structure::add_residue)
        .def("add_terminal", &model::Structure::add_terminal)
        .def("add_helix", &model::Structure::add_helix)
        .def("add_sheet", &model::Structure::add_sheet)
        .def("add_ssbond", &model::Structure::add_ssbond)
        .def("set_box_dimensions", &model::Structure::set_box_dimensions)
        .def("clear", &model::Structure::clear);
    
    // Bind Atom class
    py::class_<model::Atom, std::shared_ptr<model::Atom>>(model, "PDBAtom")
        .def(py::init<>())
        .def("get_bynu", &model::Atom::get_bynu)
        .def("get_type", &model::Atom::get_type)
        .def("get_resname", &model::Atom::get_resname)
        .def("get_ires", &model::Atom::get_ires)
        .def("get_chain", &model::Atom::get_chain)
        .def("get_coor", &model::Atom::get_coor)
        .def("is_hetatm", &model::Atom::is_hetatm)
        .def("get_occupancy", &model::Atom::get_occupancy)
        .def("get_tempfactor", &model::Atom::get_tempfactor)
        .def("get_segid", &model::Atom::get_segid)
        .def("get_inscode", &model::Atom::get_inscode)
        .def("get_formatted_atom_name", &model::Atom::get_formatted_atom_name)
        .def("get_residue_id", &model::Atom::get_residue_id)
        .def("get_mass", &model::Atom::get_mass)
        .def("set_residue_id", &model::Atom::set_residue_id)
        .def("set_coor", &model::Atom::set_coor)
        .def("set_mass_charge", &model::Atom::set_mass_charge)
        .def("set_lj_params", &model::Atom::set_lj_params)
        .def("set_occupancy", &model::Atom::set_occupancy)
        .def("set_tempfactor", &model::Atom::set_tempfactor)
        .def("set_element", &model::Atom::set_element)
        .def("set_charge_string", &model::Atom::set_charge_string)
        .def("set_chain", &model::Atom::set_chain)
        .def("set_hetatm", &model::Atom::set_hetatm)
        .def("set_bynu", &model::Atom::set_bynu)
        .def("set_type", &model::Atom::set_type)
        .def("set_resname", &model::Atom::set_resname)
        .def("set_ires", &model::Atom::set_ires)
        .def("set_segid", &model::Atom::set_segid)
        .def("set_altloc", &model::Atom::set_altloc)
        .def("set_inscode", &model::Atom::set_inscode)
        .def("has_lj_params", &model::Atom::has_lj_params)
        .def("is_valid", &model::Atom::is_valid);

    // Bind PDB Residue class
    py::class_<model::Residue, std::shared_ptr<model::Residue>>(model, "PDBResidue")
        .def(py::init<>())
        .def("get_resname", &model::Residue::get_resname)
        .def("get_ires", &model::Residue::get_ires)
        .def("get_segid", &model::Residue::get_segid)
        .def("get_iseg", &model::Residue::get_iseg)
        .def("get_chain", &model::Residue::get_chain)
        .def("get_inscode", &model::Residue::get_inscode)
        .def("get_atoms", &model::Residue::get_atoms)
        .def("add_atom", py::overload_cast<const model::Atom&>(&model::Residue::add_atom))
        .def("add_atom", py::overload_cast<std::shared_ptr<model::Atom>>(&model::Residue::add_atom))
        .def("find_atom", &model::Residue::find_atom)
        .def("atom_count", &model::Residue::atom_count)
        .def("is_valid", &model::Residue::is_valid)
        .def("calculate_center_of_mass", &model::Residue::calculate_center_of_mass)
        .def("get_center_of_mass", &model::Residue::get_center_of_mass)
        .def("has_atom_type", &model::Residue::has_atom_type)
        .def("get_residue_id", &model::Residue::get_residue_id)
        .def("set_residue_id", &model::Residue::set_residue_id)
        .def("find_atom_by_pdb_name", &model::Residue::find_atom_by_pdb_name)
        .def("update_atom_map", &model::Residue::update_atom_map)
        .def("get_atom_range", &model::Residue::get_atom_range);

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