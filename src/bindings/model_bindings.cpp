// src/bindings/model_bindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "model/atom.hpp"
#include "model/residue.hpp"
#include "model/topology.hpp"
#include "model/structure.hpp"
#include "model/forcefield.hpp"
#include "model/param.hpp"
#include "model/molecular.hpp"
#include "model/montecarlo.hpp"

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

    // Bind TopologyAtom
    py::class_<model::TopologyAtom>(model, "TopologyAtom")
        .def(py::init<>())
        .def_readwrite("id", &model::TopologyAtom::id)
        .def_readwrite("name", &model::TopologyAtom::name)
        .def_readwrite("type", &model::TopologyAtom::type)
        .def_readwrite("charge", &model::TopologyAtom::charge)
        .def_readwrite("mass", &model::TopologyAtom::mass)
        .def_readwrite("residue_id", &model::TopologyAtom::residue_id)
        .def_readwrite("segment_id", &model::TopologyAtom::segment_id)
        .def_readwrite("typeB", &model::TopologyAtom::typeB)
        .def_readwrite("chargeB", &model::TopologyAtom::chargeB)
        .def_readwrite("massB", &model::TopologyAtom::massB)
        .def_readwrite("has_b_state", &model::TopologyAtom::has_b_state);

    // Bind TopologySegment
    py::class_<model::TopologySegment>(model, "TopologySegment")
        .def(py::init<>())
        .def_readwrite("id", &model::TopologySegment::id)
        .def_readwrite("name", &model::TopologySegment::name)
        .def_readwrite("residues", &model::TopologySegment::residues);

    // Bind TopologyBond
    py::class_<model::TopologyBond>(model, "TopologyBond")
        .def(py::init<>())
        .def_readwrite("atom1", &model::TopologyBond::atom1)
        .def_readwrite("atom2", &model::TopologyBond::atom2)
        .def_readwrite("length", &model::TopologyBond::length)
        .def_readwrite("force_constant", &model::TopologyBond::force_constant)
        .def_readwrite("function_type", &model::TopologyBond::function_type)
        .def("__len__", [](const model::TopologyBond&) { return 2; })  // Bond always connects 2 atoms
        .def("__getitem__", [](const model::TopologyBond& bond, size_t i) {
            if (i == 0) return bond.atom1;
            if (i == 1) return bond.atom2;
            throw py::index_error("Bond index out of range");
        });

    // Bind TopologyAngle
    py::class_<model::TopologyAngle>(model, "TopologyAngle")
        .def(py::init<>())
        .def_readwrite("atom1", &model::TopologyAngle::atom1)
        .def_readwrite("atom2", &model::TopologyAngle::atom2)
        .def_readwrite("atom3", &model::TopologyAngle::atom3)
        .def_readwrite("angle", &model::TopologyAngle::angle)
        .def_readwrite("force_constant", &model::TopologyAngle::force_constant)
        .def_readwrite("function_type", &model::TopologyAngle::function_type)
        .def_readwrite("ub_length", &model::TopologyAngle::ub_length)
        .def_readwrite("ub_constant", &model::TopologyAngle::ub_constant)
        .def_readwrite("has_ub", &model::TopologyAngle::has_ub)
        .def("__len__", [](const model::TopologyAngle&) { return 3; })  // Angle always involves 3 atoms
        .def("__getitem__", [](const model::TopologyAngle& angle, size_t i) {
            if (i == 0) return angle.atom1;
            if (i == 1) return angle.atom2;
            if (i == 2) return angle.atom3;
            throw py::index_error("Angle index out of range");
        });

    // Bind TopologyDihedral
    py::class_<model::TopologyDihedral>(model, "TopologyDihedral")
        .def(py::init<>())
        .def_readwrite("atom1", &model::TopologyDihedral::atom1)
        .def_readwrite("atom2", &model::TopologyDihedral::atom2)
        .def_readwrite("atom3", &model::TopologyDihedral::atom3)
        .def_readwrite("atom4", &model::TopologyDihedral::atom4)
        .def_readwrite("multiplicity", &model::TopologyDihedral::multiplicity)
        .def_readwrite("angle", &model::TopologyDihedral::angle)
        .def_readwrite("force_constant", &model::TopologyDihedral::force_constant)
        .def_readwrite("improper", &model::TopologyDihedral::improper)
        .def_readwrite("function_type", &model::TopologyDihedral::function_type);

    // Bind TopologyDonor
    py::class_<model::TopologyDonor>(model, "TopologyDonor")
        .def(py::init<>())
        .def_readwrite("donor_atom", &model::TopologyDonor::donor_atom)
        .def_readwrite("hydrogen_atom", &model::TopologyDonor::hydrogen_atom);

    // Bind TopologyAcceptor
    py::class_<model::TopologyAcceptor>(model, "TopologyAcceptor")
        .def(py::init<>())
        .def_readwrite("acceptor_atom", &model::TopologyAcceptor::acceptor_atom);

    // Bind TopologyCmap
    py::class_<model::TopologyCmap>(model, "TopologyCmap")
        .def(py::init<>())
        .def_readwrite("atoms", &model::TopologyCmap::atoms)
        .def_readwrite("function_type", &model::TopologyCmap::function_type);

    // Bind StandardCmap
    py::class_<model::StandardCmap>(model, "StandardCmap")
        .def(py::init<>())
        .def_readwrite("atoms", &model::StandardCmap::atoms)
        .def_readwrite("raw_atoms", &model::StandardCmap::raw_atoms)
        .def_readwrite("is_psf_format", &model::StandardCmap::is_psf_format)
        .def_readwrite("function_type", &model::StandardCmap::function_type);

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
        .def("get_x", &model::Atom::get_x)
        .def("get_y", &model::Atom::get_y)
        .def("get_z", &model::Atom::get_z)
        .def("is_hetatm", &model::Atom::is_hetatm)
        .def("get_occupancy", &model::Atom::get_occupancy)
        .def("get_tempfactor", &model::Atom::get_tempfactor)
        .def("get_segid", &model::Atom::get_segid)
        .def("get_inscode", &model::Atom::get_inscode)
        .def("get_formatted_atom_name", &model::Atom::get_formatted_atom_name)
        .def("get_residue_id", &model::Atom::get_residue_id)
        .def("get_mass", &model::Atom::get_mass)
        .def("get_charge", &model::Atom::get_charge)
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
        .def("set_resname", &model::Residue::set_resname)
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
        .def("has_atom", &model::Topology::has_atom)
        .def("has_residue", &model::Topology::has_residue)
        .def("has_segment", &model::Topology::has_segment)
        .def("has_bond", &model::Topology::has_bond)
        .def("has_angle", &model::Topology::has_angle)
        .def("has_dihedral", &model::Topology::has_dihedral)
        .def("has_improper", &model::Topology::has_improper)
        .def("has_donor", &model::Topology::has_donor)
        .def("has_acceptor", &model::Topology::has_acceptor)
        .def("has_cmap", &model::Topology::has_cmap)
        .def("has_group", &model::Topology::has_group)
        .def("get_atom", &model::Topology::get_atom, py::return_value_policy::reference_internal, "Get atom by index")
        .def("get_residue", &model::Topology::get_residue, py::return_value_policy::reference_internal, "Get residue by index")
        .def("get_segment", &model::Topology::get_segment, py::return_value_policy::reference_internal, "Get segment by index")
        .def("get_group", &model::Topology::get_group, py::return_value_policy::reference_internal, "Get group by index")
        .def("find_atom", &model::Topology::find_atom)
        .def("find_residue", &model::Topology::find_residue)
        .def("find_segment", &model::Topology::find_segment)
        .def("get_cmaps", &model::Topology::get_cmaps, py::return_value_policy::reference_internal, "Get all CMAP terms")
        .def_property_readonly("cmaps", [](const model::Topology& t) { return t.get_cmaps(); }, "Get all CMAP terms");

    // Bind TopologyResidue
    py::class_<model::TopologyResidue>(model, "TopologyResidue")
        .def_readonly("name", &model::TopologyResidue::name)
        .def_readonly("number", &model::TopologyResidue::number)
        .def_readonly("atoms", &model::TopologyResidue::atoms)
        .def_readonly("segment", &model::TopologyResidue::segment);

    // Bind TopologyGroup
    py::class_<model::TopologyGroup>(model, "TopologyGroup")
        .def_readonly("id", &model::TopologyGroup::id)
        .def_readonly("atoms", &model::TopologyGroup::atoms)
        .def_readonly("type", &model::TopologyGroup::type);

    // NonbondedParams
    py::class_<model::NonbondedParams>(m, "NonbondedParams")
        .def(py::init<>())
        .def_readwrite("nbxmod", &model::NonbondedParams::nbxmod)
        .def_readwrite("cdiel", &model::NonbondedParams::cdiel)
        .def_readwrite("fshift", &model::NonbondedParams::fshift)
        .def_readwrite("vatom", &model::NonbondedParams::vatom)
        .def_readwrite("vdistance", &model::NonbondedParams::vdistance)
        .def_readwrite("vfswitch", &model::NonbondedParams::vfswitch)
        .def_readwrite("cutnb", &model::NonbondedParams::cutnb)
        .def_readwrite("ctofnb", &model::NonbondedParams::ctofnb)
        .def_readwrite("ctonnb", &model::NonbondedParams::ctonnb)
        .def_readwrite("eps", &model::NonbondedParams::eps)
        .def_readwrite("e14fac", &model::NonbondedParams::e14fac)
        .def_readwrite("wmin", &model::NonbondedParams::wmin);

    // LJParams
    py::class_<model::LJParams>(m, "LJParams")
        .def(py::init<>())
        .def_readwrite("epsilon", &model::LJParams::epsilon)
        .def_readwrite("rmin_half", &model::LJParams::rmin_half);

    // BondParams
    py::class_<model::BondParams>(m, "BondParams")
        .def(py::init<>())
        .def_readwrite("kb", &model::BondParams::kb)
        .def_readwrite("b0", &model::BondParams::b0);

    // AngleParams
    py::class_<model::AngleParams>(m, "AngleParams")
        .def(py::init<>())
        .def_readwrite("ktheta", &model::AngleParams::ktheta)
        .def_readwrite("theta0", &model::AngleParams::theta0)
        .def_readwrite("kub", &model::AngleParams::kub)
        .def_readwrite("s0", &model::AngleParams::s0);

    // DihedralParams
    py::class_<model::DihedralParams>(m, "DihedralParams")
        .def(py::init<>())
        .def_readwrite("kchi", &model::DihedralParams::kchi)
        .def_readwrite("n", &model::DihedralParams::n)
        .def_readwrite("delta", &model::DihedralParams::delta);

    // ImproperParams
    py::class_<model::ImproperParams>(m, "ImproperParams")
        .def(py::init<>())
        .def_readwrite("kpsi", &model::ImproperParams::kpsi)
        .def_readwrite("psi0", &model::ImproperParams::psi0);

    // NBFIXParams
    py::class_<model::NBFIXParams>(m, "NBFIXParams")
        .def(py::init<>())
        .def_readwrite("epsilon", &model::NBFIXParams::epsilon)
        .def_readwrite("rmin", &model::NBFIXParams::rmin);

    // Bind ForceField
    py::class_<model::ForceField>(m, "ForceField")
        .def(py::init<>())
        .def("add_atom_mass", &model::ForceField::add_atom_mass)
        .def("add_lj_params", &model::ForceField::add_lj_params)
        .def("add_nbfix", &model::ForceField::add_nbfix)
        .def("add_bond_params", &model::ForceField::add_bond_params)
        .def("add_angle_params", &model::ForceField::add_angle_params)
        .def("add_dihedral_params", &model::ForceField::add_dihedral_params)
        .def("add_improper_params", &model::ForceField::add_improper_params)
        .def("get_atom_mass", &model::ForceField::get_atom_mass)
        .def("get_lj_params", static_cast<const model::LJParams& (model::ForceField::*)(const std::string&) const>(&model::ForceField::get_lj_params))
        .def("get_nbfix", [](const model::ForceField& ff, const std::string& type1, const std::string& type2) {
            auto result = ff.get_nbfix(type1, type2);
            return std::make_tuple(result.first.epsilon, result.first.rmin, result.second);
        }, "Get NBFIX parameters for a pair of atom types. Returns (epsilon, rmin, found)")
        .def("get_bond_params", static_cast<const model::BondParams& (model::ForceField::*)(const std::string&, const std::string&) const>(&model::ForceField::get_bond_params))
        .def("get_angle_params", static_cast<const model::AngleParams& (model::ForceField::*)(const std::string&, const std::string&, const std::string&) const>(&model::ForceField::get_angle_params))
        .def("get_dihedral_params", static_cast<const std::vector<model::DihedralParams>& (model::ForceField::*)(const std::string&, const std::string&, const std::string&, const std::string&) const>(&model::ForceField::get_dihedral_params))
        .def("get_improper_params", static_cast<const model::ImproperParams& (model::ForceField::*)(const std::string&, const std::string&, const std::string&, const std::string&) const>(&model::ForceField::get_improper_params))
        .def("has_atom_mass", &model::ForceField::has_atom_mass)
        .def("has_lj_params", &model::ForceField::has_lj_params)
        .def("has_nbfix", &model::ForceField::has_nbfix)
        .def("has_bond_params", &model::ForceField::has_bond_params)
        .def("has_angle_params", &model::ForceField::has_angle_params)
        .def("has_dihedral_params", &model::ForceField::has_dihedral_params)
        .def("has_improper_params", &model::ForceField::has_improper_params)
        .def("get_num_atom_types", &model::ForceField::get_num_atom_types)
        .def("get_num_lj_params", &model::ForceField::get_num_lj_params)
        .def("get_num_nbfix", &model::ForceField::get_num_nbfix)
        .def("get_num_bond_types", &model::ForceField::get_num_bond_types)
        .def("get_num_angle_types", &model::ForceField::get_num_angle_types)
        .def("get_num_dihedral_types", &model::ForceField::get_num_dihedral_types)
        .def("get_num_improper_types", &model::ForceField::get_num_improper_types)
        .def("get_nonbonded_params", static_cast<const model::NonbondedParams& (model::ForceField::*)() const>(&model::ForceField::get_nonbonded_params))
        .def_static("makeTypePair", &model::ForceField::makeTypePair)
        .def_static("makeTypeTriple", &model::ForceField::makeTypeTriple)
        .def_static("makeTypeQuad", &model::ForceField::makeTypeQuad)
        // Property accessors
        .def_property_readonly("atom_masses", static_cast<const std::map<std::string, double>& (model::ForceField::*)() const>(&model::ForceField::get_atom_masses))
        .def_property_readonly("lj_params", static_cast<const std::map<std::string, model::LJParams>& (model::ForceField::*)() const>(&model::ForceField::get_lj_params))
        .def_property_readonly("nbfix", static_cast<const std::map<std::pair<std::string, std::string>, model::NBFIXParams>& (model::ForceField::*)() const>(&model::ForceField::get_nbfix))
        .def_property_readonly("bond_params", static_cast<const std::map<std::pair<std::string, std::string>, model::BondParams>& (model::ForceField::*)() const>(&model::ForceField::get_bond_params))
        .def_property_readonly("angle_params", static_cast<const std::map<std::tuple<std::string, std::string, std::string>, model::AngleParams>& (model::ForceField::*)() const>(&model::ForceField::get_angle_params))
        .def_property_readonly("dihedral_params", static_cast<const std::map<std::tuple<std::string, std::string, std::string, std::string>, std::vector<model::DihedralParams>>& (model::ForceField::*)() const>(&model::ForceField::get_dihedral_params))
        .def_property_readonly("improper_params", static_cast<const std::map<std::tuple<std::string, std::string, std::string, std::string>, model::ImproperParams>& (model::ForceField::*)() const>(&model::ForceField::get_improper_params));

    // Bind Param class and its nested structs
    auto param = py::class_<model::Param>(model, "Param")
        .def(py::init<>())
        // Basic info
        .def_property("basic_info",
            py::overload_cast<>(&model::Param::get_basic_info, py::const_),
            py::overload_cast<>(&model::Param::get_basic_info))
        // Space info
        .def_property("space_info",
            py::overload_cast<>(&model::Param::get_space_info, py::const_),
            py::overload_cast<>(&model::Param::get_space_info))
        // MC info
        .def_property("mc_info",
            py::overload_cast<>(&model::Param::get_mc_info, py::const_),
            py::overload_cast<>(&model::Param::get_mc_info))
        // Energy info
        .def_property("energy_info",
            py::overload_cast<>(&model::Param::get_energy_info, py::const_),
            py::overload_cast<>(&model::Param::get_energy_info))
        // Fragment info
        .def_property("fragment_info",
            py::overload_cast<>(&model::Param::get_fragment_info, py::const_),
            py::overload_cast<>(&model::Param::get_fragment_info))
        // Bias info
        .def_property("bias_info",
            py::overload_cast<>(&model::Param::get_bias_info, py::const_),
            py::overload_cast<>(&model::Param::get_bias_info))
        // File info
        .def_property("file_info",
            py::overload_cast<>(&model::Param::get_file_info, py::const_),
            py::overload_cast<>(&model::Param::get_file_info))
        .def("clear", &model::Param::clear);

    // Bind BasicInfo struct
    py::class_<model::Param::BasicInfo>(param, "BasicInfo")
        .def(py::init<>())
        .def_readwrite("version", &model::Param::BasicInfo::version)
        .def_readwrite("verbosity", &model::Param::BasicInfo::verbosity)
        .def_readwrite("debug", &model::Param::BasicInfo::debug)
        .def_readwrite("print_logfile", &model::Param::BasicInfo::print_logfile)
        .def_readwrite("param_file", &model::Param::BasicInfo::param_file)
        .def_readwrite("log_file", &model::Param::BasicInfo::log_file)
        .def_readwrite("random_seed", &model::Param::BasicInfo::random_seed)
        .def_readwrite("num_threads", &model::Param::BasicInfo::num_threads)
        .def_readwrite("is_box", &model::Param::BasicInfo::is_box)
        .def_readwrite("init_cycle", &model::Param::BasicInfo::init_cycle)
        .def_readwrite("conserve_fragments", &model::Param::BasicInfo::conserve_fragments);

    // Bind SpaceInfo struct
    py::class_<model::Param::SpaceInfo>(param, "SpaceInfo")
        .def(py::init<>())
        .def_readwrite("grid_spacing", &model::Param::SpaceInfo::grid_spacing)
        .def_readwrite("gc_center", &model::Param::SpaceInfo::gc_center)
        .def_readwrite("sys_center", &model::Param::SpaceInfo::sys_center)
        .def_readwrite("crystal_dim", &model::Param::SpaceInfo::crystal_dim)
        .def_readwrite("box_size", &model::Param::SpaceInfo::box_size)
        .def_readwrite("volume", &model::Param::SpaceInfo::volume)
        .def_readwrite("target_volume", &model::Param::SpaceInfo::target_volume)
        .def_readwrite("sys_box_volume", &model::Param::SpaceInfo::sys_box_volume)
        .def_readwrite("gcmc_volume", &model::Param::SpaceInfo::gcmc_volume)
        .def_readwrite("protein_volume", &model::Param::SpaceInfo::protein_volume)
        .def_readwrite("use_vdw_radius_for_grid", &model::Param::SpaceInfo::use_vdw_radius_for_grid)
        .def_readwrite("exclude_hydrogens_from_grid", &model::Param::SpaceInfo::exclude_hydrogens_from_grid)
        .def_readwrite("exclude_protein_volume", &model::Param::SpaceInfo::exclude_protein_volume)
        .def_readwrite("tmp_prob", &model::Param::SpaceInfo::tmp_prob)
        .def_readwrite("cutoff", &model::Param::SpaceInfo::cutoff);

    // Bind MCInfo struct
    py::class_<model::Param::MCInfo>(param, "MCInfo")
        .def(py::init<>())
        .def_readwrite("mc_steps", &model::Param::MCInfo::mc_steps)
        .def_readwrite("current_step", &model::Param::MCInfo::current_step)
        .def_readwrite("print_freq", &model::Param::MCInfo::print_freq)
        .def_readwrite("temperature", &model::Param::MCInfo::temperature)
        .def_readwrite("beta", &model::Param::MCInfo::beta)
        .def_readwrite("insertion_deletion_frac", &model::Param::MCInfo::insertion_deletion_frac)
        .def_readwrite("translation_rotation_frac", &model::Param::MCInfo::translation_rotation_frac)
        .def_readwrite("max_translation_dist", &model::Param::MCInfo::max_translation_dist)
        .def_readwrite("max_rotation_angle", &model::Param::MCInfo::max_rotation_angle)
        .def_readwrite("operation_types", &model::Param::MCInfo::operation_types)
        .def_readwrite("mc_time_list", &model::Param::MCInfo::mc_time_list)
        .def_readwrite("mc_time_cumulative", &model::Param::MCInfo::mc_time_cumulative)
        .def_readwrite("fragment_prob", &model::Param::MCInfo::fragment_prob)
        .def_readwrite("water_prob", &model::Param::MCInfo::water_prob)
        .def_readwrite("atom_prob", &model::Param::MCInfo::atom_prob)
        .def_readwrite("test_prob", &model::Param::MCInfo::test_prob)
        .def_readwrite("rotate_dih_status", &model::Param::MCInfo::rotate_dih_status)
        .def_readonly("BOLTZMANN", &model::Param::MCInfo::BOLTZMANN)
        .def_readonly("KCAL_TO_KJ", &model::Param::MCInfo::KCAL_TO_KJ);

    // Bind EnergyInfo struct
    py::class_<model::Param::EnergyInfo>(param, "EnergyInfo")
        .def(py::init<>())
        .def_readwrite("use_group_cutoff", &model::Param::EnergyInfo::use_group_cutoff)
        .def_readwrite("fragment_cutoff", &model::Param::EnergyInfo::fragment_cutoff)
        .def_readwrite("protein_cutoff", &model::Param::EnergyInfo::protein_cutoff)
        .def_readwrite("fragment_cutoff_squared", &model::Param::EnergyInfo::fragment_cutoff_squared)
        .def_readwrite("protein_cutoff_squared", &model::Param::EnergyInfo::protein_cutoff_squared)
        .def_readwrite("pairlist_cutoff", &model::Param::EnergyInfo::pairlist_cutoff)
        .def_readwrite("pairlist_cutoff_squared", &model::Param::EnergyInfo::pairlist_cutoff_squared)
        .def_readwrite("pairlist_freq", &model::Param::EnergyInfo::pairlist_freq)
        .def_readwrite("use_switching", &model::Param::EnergyInfo::use_switching)
        .def_readwrite("switch_dist_fragment", &model::Param::EnergyInfo::switch_dist_fragment)
        .def_readwrite("switch_dist_protein", &model::Param::EnergyInfo::switch_dist_protein)
        .def_readwrite("switch_dist_fragment_squared", &model::Param::EnergyInfo::switch_dist_fragment_squared)
        .def_readwrite("switch_dist_protein_squared", &model::Param::EnergyInfo::switch_dist_protein_squared)
        .def_readwrite("energy_sw_ref", &model::Param::EnergyInfo::energy_sw_ref)
        .def_readwrite("energy_sw_scale", &model::Param::EnergyInfo::energy_sw_scale)
        .def_readwrite("test_sw_filters", &model::Param::EnergyInfo::test_sw_filters)
        .def_readwrite("apply_sw_filters", &model::Param::EnergyInfo::apply_sw_filters)
        .def_readwrite("test_energy", &model::Param::EnergyInfo::test_energy)
        .def_readwrite("pair_list_cutoff_fragment", &model::Param::EnergyInfo::pair_list_cutoff_fragment)
        .def_readwrite("pair_list_cutoff_protein", &model::Param::EnergyInfo::pair_list_cutoff_protein)
        .def_readwrite("pair_list_cutoff_fragment_squared", &model::Param::EnergyInfo::pair_list_cutoff_fragment_squared)
        .def_readwrite("pair_list_cutoff_protein_squared", &model::Param::EnergyInfo::pair_list_cutoff_protein_squared);

    // Bind FragmentInfo struct
    py::class_<model::Param::FragmentInfo>(param, "FragmentInfo")
        .def(py::init<>())
        .def_readwrite("water_density", &model::Param::FragmentInfo::water_density)
        .def_readwrite("epsilon", &model::Param::FragmentInfo::epsilon)
        .def_readwrite("num_waters", &model::Param::FragmentInfo::num_waters)
        .def_readwrite("target_num_waters", &model::Param::FragmentInfo::target_num_waters)
        .def_readwrite("water_index", &model::Param::FragmentInfo::water_index)
        .def_readwrite("excess_threshold", &model::Param::FragmentInfo::excess_threshold)
        .def_readwrite("use_number_water_nbar", &model::Param::FragmentInfo::use_number_water_nbar)
        .def_readwrite("use_const_water_nbar", &model::Param::FragmentInfo::use_const_water_nbar)
        .def_readwrite("const_water_nbar", &model::Param::FragmentInfo::const_water_nbar)
        .def_readwrite("init_cutoff", &model::Param::FragmentInfo::init_cutoff)
        .def_readwrite("init_cutoff_squared", &model::Param::FragmentInfo::init_cutoff_squared)
        .def_readwrite("use_gcmc_cutoff", &model::Param::FragmentInfo::use_gcmc_cutoff)
        .def_readwrite("gcmc_cutoff", &model::Param::FragmentInfo::gcmc_cutoff)
        .def_readwrite("gcmc_cutoff_squared", &model::Param::FragmentInfo::gcmc_cutoff_squared)
        .def_readwrite("remove_init", &model::Param::FragmentInfo::remove_init)
        .def_readwrite("remove_excess", &model::Param::FragmentInfo::remove_excess)
        .def_readwrite("confs_list", &model::Param::FragmentInfo::confs_list)
        .def_readwrite("cavity_index_list", &model::Param::FragmentInfo::cavity_index_list)
        .def_readwrite("cavity_list", &model::Param::FragmentInfo::cavity_list)
        .def_readwrite("conc_list", &model::Param::FragmentInfo::conc_list)
        .def_readwrite("muex_list", &model::Param::FragmentInfo::muex_list)
        .def_readwrite("radius_list", &model::Param::FragmentInfo::radius_list)
        .def_readwrite("conf_list", &model::Param::FragmentInfo::conf_list)
        .def_readwrite("flag_remove_init", &model::Param::FragmentInfo::flag_remove_init)
        .def_readwrite("flag_remove_excess", &model::Param::FragmentInfo::flag_remove_excess)
        .def_readwrite("total_protitp_size", &model::Param::FragmentInfo::total_protitp_size)
        .def_readwrite("fragconf_list", &model::Param::FragmentInfo::fragconf_list);

    // Bind BiasInfo struct
    py::class_<model::Param::BiasInfo>(param, "BiasInfo")
        .def(py::init<>())
        .def_readwrite("use_cavity_bias", &model::Param::BiasInfo::use_cavity_bias)
        .def_readwrite("sigma", &model::Param::BiasInfo::sigma)
        .def_readwrite("sigma_squared", &model::Param::BiasInfo::sigma_squared)
        .def_readwrite("use_conf_bias", &model::Param::BiasInfo::use_conf_bias)
        .def_readwrite("num_conf_bias_trials", &model::Param::BiasInfo::num_conf_bias_trials);

    // Bind FileInfo struct
    py::class_<model::Param::FileInfo>(param, "FileInfo")
        .def(py::init<>())
        .def_readwrite("topology_file", &model::Param::FileInfo::topology_file)
        .def_readwrite("input_pdb_file", &model::Param::FileInfo::input_pdb_file)
        .def_readwrite("output_pdb_file", &model::Param::FileInfo::output_pdb_file)
        .def_readwrite("output_top_file", &model::Param::FileInfo::output_top_file)
        .def_readwrite("atomtype_file", &model::Param::FileInfo::atomtype_file)
        .def_readwrite("monomer_dir", &model::Param::FileInfo::monomer_dir)
        .def_readwrite("conc_norm", &model::Param::FileInfo::conc_norm)
        .def_readwrite("conc_region", &model::Param::FileInfo::conc_region)
        .def_readwrite("par_files", &model::Param::FileInfo::par_files)
        .def_readwrite("protein_top_files", &model::Param::FileInfo::protein_top_files)
        .def_readwrite("fragment_top_files", &model::Param::FileInfo::fragment_top_files)
        .def_readwrite("fragment_names", &model::Param::FileInfo::fragment_names)
        .def_readwrite("fragment_mqtr_files", &model::Param::FileInfo::fragment_mqtr_files)
        .def_readwrite("tmp_frag_name", &model::Param::FileInfo::tmp_frag_name)
        .def_readwrite("generate_maps", &model::Param::FileInfo::generate_maps)
        .def_readwrite("map_prefix", &model::Param::FileInfo::map_prefix);

    // Bind Molecular class
    py::class_<model::Molecular, std::shared_ptr<model::Molecular>>(model, "Molecular")
        .def(py::init<>())
        .def_readwrite("atoms", &model::Molecular::atoms)
        .def_readwrite("residues", &model::Molecular::residues)
        .def_readwrite("terminals", &model::Molecular::terminals)
        .def_readwrite("helices", &model::Molecular::helices)
        .def_readwrite("sheets", &model::Molecular::sheets)
        .def_readwrite("ssbonds", &model::Molecular::ssbonds)
        .def_readwrite("boxDimensions", &model::Molecular::boxDimensions)
        .def_readwrite("topology_atoms", &model::Molecular::topology_atoms)
        .def_readwrite("topology_residues", &model::Molecular::topology_residues)
        .def_readwrite("segments", &model::Molecular::segments)
        .def_readwrite("bonds", &model::Molecular::bonds)
        .def_readwrite("angles", &model::Molecular::angles)
        .def_readwrite("dihedrals", &model::Molecular::dihedrals)
        .def_readwrite("donors", &model::Molecular::donors)
        .def_readwrite("acceptors", &model::Molecular::acceptors)
        .def_readwrite("exclusions", &model::Molecular::exclusions)
        .def_readwrite("groups", &model::Molecular::groups)
        .def_readwrite("cmaps", &model::Molecular::cmaps)
        .def_readwrite("standard_cmaps", &model::Molecular::standard_cmaps)
        .def_readwrite("titles", &model::Molecular::titles)
        .def_readwrite("segment_map", &model::Molecular::segment_map)
        .def_readwrite("residue_map", &model::Molecular::residue_map)
        .def_readwrite("atom_map", &model::Molecular::atom_map)
        .def("get_num_atoms", &model::Molecular::get_num_atoms)
        .def("get_num_residues", &model::Molecular::get_num_residues)
        .def("get_num_segments", &model::Molecular::get_num_segments)
        .def("get_num_bonds", &model::Molecular::get_num_bonds)
        .def("get_num_angles", &model::Molecular::get_num_angles)
        .def("get_num_dihedrals", &model::Molecular::get_num_dihedrals)
        .def("get_num_impropers", &model::Molecular::get_num_impropers)
        .def("get_num_standard_cmaps", &model::Molecular::get_num_standard_cmaps)
        .def("add_standard_cmap", &model::Molecular::add_standard_cmap)
        .def("clear", &model::Molecular::clear);

    // Bind GCMCInfo
    py::class_<pygcmc::model::MCInfo>(m, "MCInfo")
        .def(py::init<>())
        .def_readwrite("mc_steps", &pygcmc::model::MCInfo::mcSteps)
        .def_property("box",
            [](const pygcmc::model::MCInfo& info) {
                return std::vector<float>{info.box[0], info.box[1], info.box[2]};
            },
            [](pygcmc::model::MCInfo& info, const std::vector<float>& box) {
                if (box.size() != 3) throw std::runtime_error("Box must have 3 dimensions");
                info.box[0] = box[0];
                info.box[1] = box[1];
                info.box[2] = box[2];
            })
        .def_readwrite("cutoff", &pygcmc::model::MCInfo::cutoff)
        .def_readwrite("beta", &pygcmc::model::MCInfo::beta)
        .def_readwrite("max_residues", &pygcmc::model::MCInfo::maxResidues)
        .def_readwrite("max_atoms", &pygcmc::model::MCInfo::maxAtoms)
        .def_readwrite("max_types", &pygcmc::model::MCInfo::maxTypes)
        .def_readwrite("volume", &pygcmc::model::MCInfo::volume)
        .def_readwrite("seed", &pygcmc::model::MCInfo::seed)
        .def("setTemperature", &pygcmc::model::MCInfo::setTemperature, "Set temperature in Kelvin and calculate beta");

    // Bind GCMCInfo::Statistics
    py::class_<pygcmc::model::MCInfo::Statistics>(m, "MCStatistics")
        .def(py::init<>())
        .def_readwrite("totalMoves", &pygcmc::model::MCInfo::Statistics::totalMoves)
        .def_readwrite("acceptedMoves", &pygcmc::model::MCInfo::Statistics::acceptedMoves)
        .def_readwrite("insertionAttempts", &pygcmc::model::MCInfo::Statistics::insertionAttempts)
        .def_readwrite("acceptedInsertions", &pygcmc::model::MCInfo::Statistics::acceptedInsertions)
        .def_readwrite("deletionAttempts", &pygcmc::model::MCInfo::Statistics::deletionAttempts)
        .def_readwrite("acceptedDeletions", &pygcmc::model::MCInfo::Statistics::acceptedDeletions);

    // Bind MCResidue
    py::class_<pygcmc::model::MCResidue>(m, "MCResidue")
        .def(py::init<>())
        .def_readwrite("atomStart", &pygcmc::model::MCResidue::atomStart)
        .def_property_readonly("atom_start", [](const pygcmc::model::MCResidue& r) { return r.atomStart; })
        .def_readwrite("atomCount", &pygcmc::model::MCResidue::atomCount)
        .def_property_readonly("atom_count", [](const pygcmc::model::MCResidue& r) { return r.atomCount; })
        .def_readwrite("active", &pygcmc::model::MCResidue::active)
        .def_readwrite("fixed", &pygcmc::model::MCResidue::fixed)
        .def_property("center",
            [](const pygcmc::model::MCResidue& res) -> std::vector<float> {
                return {res.center[0], res.center[1], res.center[2]};
            },
            [](pygcmc::model::MCResidue& res, const std::vector<float>& center) {
                if (center.size() != 3) {
                    throw std::runtime_error("Center must be a vector of 3 floats");
                }
                res.center[0] = center[0];
                res.center[1] = center[1];
                res.center[2] = center[2];
            })
        .def_readwrite("concentration", &pygcmc::model::MCResidue::concentration)
        .def_readwrite("chemPot", &pygcmc::model::MCResidue::chemPot)
        .def_property_readonly("chem_pot", [](const pygcmc::model::MCResidue& r) { return r.chemPot; })
        .def_readwrite("type", &pygcmc::model::MCResidue::type)
        .def_readwrite("radius", &pygcmc::model::MCResidue::radius)
        .def_readwrite("energy_vdw", &pygcmc::model::MCResidue::energy_vdw)
        .def_readwrite("energy_elec", &pygcmc::model::MCResidue::energy_elec);

    // Bind MCAtom
    py::class_<pygcmc::model::MCAtom>(m, "MCAtom")
        .def(py::init<>())
        .def_readwrite("x", &pygcmc::model::MCAtom::x)
        .def_readwrite("y", &pygcmc::model::MCAtom::y)
        .def_readwrite("z", &pygcmc::model::MCAtom::z)
        .def_readwrite("charge", &pygcmc::model::MCAtom::charge)
        .def_readwrite("type", &pygcmc::model::MCAtom::type);

    // Bind MCMovementResidueInfo
    py::class_<pygcmc::model::MCMovementResidueInfo>(m, "MCMovementResidueInfo")
        .def(py::init<>())
        .def_readwrite("startIndex", &pygcmc::model::MCMovementResidueInfo::startIndex)
        .def_readwrite("activeCount", &pygcmc::model::MCMovementResidueInfo::activeCount)
        .def_readwrite("totalCount", &pygcmc::model::MCMovementResidueInfo::totalCount)
        .def_readwrite("resName", &pygcmc::model::MCMovementResidueInfo::resName);

    // Bind MCState
    py::class_<pygcmc::model::MCState>(m, "MCState")
        .def(py::init<>())
        .def("copy", [](const pygcmc::model::MCState& state) {
            pygcmc::model::MCState new_state;
            new_state.atoms = state.atoms;
            new_state.residues = state.residues;
            new_state.residueTypes = state.residueTypes;
            new_state.atomTypes = state.atomTypes;
            new_state.activeAtomCount = state.activeAtomCount;
            new_state.activeResidueCount = state.activeResidueCount;
            new_state.info = state.info;
            new_state.forcefield = state.forcefield;
            new_state.movementResidues = state.movementResidues;
            new_state.movementAtomTypes = state.movementAtomTypes;
            new_state.numMovementAtomTypes = state.numMovementAtomTypes;
            return new_state;
        }, "Create a deep copy of the MCState object")
        .def_property("atoms",
            [](pygcmc::model::MCState& state) -> std::vector<std::reference_wrapper<pygcmc::model::MCAtom>> {
                std::vector<std::reference_wrapper<pygcmc::model::MCAtom>> refs;
                refs.reserve(state.activeAtomCount);
                for (int i = 0; i < state.activeAtomCount; ++i) {
                    refs.push_back(std::ref(state.atoms[i]));
                }
                return refs;
            },
            [](pygcmc::model::MCState& state, const std::vector<pygcmc::model::MCAtom>& atoms) {
                state.atoms = atoms;
            })
        .def_property("residues",
            [](const pygcmc::model::MCState& state) {
                return std::vector<pygcmc::model::MCResidue>(state.residues.begin(), 
                    state.residues.begin() + state.activeResidueCount);
            },
            [](pygcmc::model::MCState& state, const std::vector<pygcmc::model::MCResidue>& residues) {
                state.residues = residues;
            })
        .def_readwrite("residueTypes", &pygcmc::model::MCState::residueTypes)
        .def_readwrite("atomTypes", &pygcmc::model::MCState::atomTypes)
        .def_readwrite("activeAtomCount", &pygcmc::model::MCState::activeAtomCount)
        .def_readwrite("activeResidueCount", &pygcmc::model::MCState::activeResidueCount)
        .def_readwrite("info", &pygcmc::model::MCState::info)
        .def_readwrite("forcefield", &pygcmc::model::MCState::forcefield)
        .def_property("movementResidues",
            [](const pygcmc::model::MCState& state) {
                return state.movementResidues;
            },
            [](pygcmc::model::MCState& state, const std::vector<pygcmc::model::MCMovementResidueInfo>& movementResidues) {
                state.movementResidues = movementResidues;
            })
        .def_readwrite("movementAtomTypes", &pygcmc::model::MCState::movementAtomTypes)
        .def_readwrite("numMovementAtomTypes", &pygcmc::model::MCState::numMovementAtomTypes);

    // Bind MCForceField
    py::class_<pygcmc::model::MCForceField>(m, "MCForceField")
        .def(py::init<>())
        .def_readwrite("numTotalTypes", &pygcmc::model::MCForceField::numTotalTypes)
        .def_property("maxTypes",
            [](const pygcmc::model::MCForceField& ff) { return ff.numTotalTypes; },
            [](pygcmc::model::MCForceField& ff, int value) { ff.numTotalTypes = value; })
        .def_readwrite("numMovementTypes", &pygcmc::model::MCForceField::numMovementTypes)
        .def_readwrite("ljSigma", &pygcmc::model::MCForceField::ljSigma)
        .def_readwrite("ljEps", &pygcmc::model::MCForceField::ljEps);

    // Add COULOMB constant to the module
    m.attr("COULOMB") = 138.935458; // kJ·mol^-1·nm·e^-2, Coulomb's constant in MD units
}

} // namespace bindings
} // namespace pygcmc 