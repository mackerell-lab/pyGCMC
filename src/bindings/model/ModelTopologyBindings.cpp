// src/bindings/model/ModelTopologyBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "model/ModelModule.hpp"

namespace py = pybind11;
using namespace pygcmc::model;
using namespace pygcmc::model::topology;

namespace pygcmc {
namespace bindings {
namespace model {

void init_topology_bindings(py::module&, py::module& model_module) {
    // Bind TopologyAtom
    py::class_<::pygcmc::model::TopologyAtom>(model_module, "TopologyAtom")
        .def(py::init<>())
        .def_readwrite("id", &::pygcmc::model::TopologyAtom::id)
        .def_readwrite("name", &::pygcmc::model::TopologyAtom::name)
        .def_readwrite("type", &::pygcmc::model::TopologyAtom::type)
        .def_readwrite("charge", &::pygcmc::model::TopologyAtom::charge)
        .def_readwrite("mass", &::pygcmc::model::TopologyAtom::mass)
        .def_readwrite("residue_id", &::pygcmc::model::TopologyAtom::residue_id)
        .def_readwrite("segment_id", &::pygcmc::model::TopologyAtom::segment_id)
        .def_readwrite("alpha", &::pygcmc::model::TopologyAtom::alpha)
        .def_readwrite("thole", &::pygcmc::model::TopologyAtom::thole)
        .def("get_alpha", &::pygcmc::model::TopologyAtom::get_alpha)
        .def("get_thole", &::pygcmc::model::TopologyAtom::get_thole)
        .def("set_alpha", &::pygcmc::model::TopologyAtom::set_alpha)
        .def("set_thole", &::pygcmc::model::TopologyAtom::set_thole)
        .def("set_drude_params", &::pygcmc::model::TopologyAtom::set_drude_params);

    // Bind TopologySegment
    py::class_<::pygcmc::model::TopologySegment>(model_module, "TopologySegment")
        .def(py::init<>())
        .def_readwrite("id", &::pygcmc::model::TopologySegment::id)
        .def_readwrite("name", &::pygcmc::model::TopologySegment::name)
        .def_readwrite("residues", &::pygcmc::model::TopologySegment::residues);

    // Bind TopologyBond
    py::class_<::pygcmc::model::TopologyBond>(model_module, "TopologyBond")
        .def(py::init<>())
        .def_readwrite("atom1", &::pygcmc::model::TopologyBond::atom1)
        .def_readwrite("atom2", &::pygcmc::model::TopologyBond::atom2)
        .def_readwrite("length", &::pygcmc::model::TopologyBond::length)
        .def_readwrite("force_constant", &::pygcmc::model::TopologyBond::force_constant)
        .def_readwrite("function_type", &::pygcmc::model::TopologyBond::function_type)
        .def("__len__", [](const ::pygcmc::model::TopologyBond&) { return 2; })  // Bond always connects 2 atoms
        .def("__getitem__", [](const ::pygcmc::model::TopologyBond& bond, size_t i) {
            if (i == 0) return bond.atom1;
            if (i == 1) return bond.atom2;
            throw py::index_error("Bond index out of range");
        });

    // Bind TopologyAngle
    py::class_<::pygcmc::model::TopologyAngle>(model_module, "TopologyAngle")
        .def(py::init<>())
        .def_readwrite("atom1", &::pygcmc::model::TopologyAngle::atom1)
        .def_readwrite("atom2", &::pygcmc::model::TopologyAngle::atom2)
        .def_readwrite("atom3", &::pygcmc::model::TopologyAngle::atom3)
        .def_readwrite("angle", &::pygcmc::model::TopologyAngle::angle)
        .def_readwrite("force_constant", &::pygcmc::model::TopologyAngle::force_constant)
        .def_readwrite("function_type", &::pygcmc::model::TopologyAngle::function_type)
        .def("__len__", [](const ::pygcmc::model::TopologyAngle&) { return 3; })  // Angle always involves 3 atoms
        .def("__getitem__", [](const ::pygcmc::model::TopologyAngle& angle, size_t i) {
            if (i == 0) return angle.atom1;
            if (i == 1) return angle.atom2;
            if (i == 2) return angle.atom3;
            throw py::index_error("Angle index out of range");
        });

    // Bind TopologyDihedral
    py::class_<::pygcmc::model::TopologyDihedral>(model_module, "TopologyDihedral")
        .def(py::init<>())
        .def_readwrite("atom1", &::pygcmc::model::TopologyDihedral::atom1)
        .def_readwrite("atom2", &::pygcmc::model::TopologyDihedral::atom2)
        .def_readwrite("atom3", &::pygcmc::model::TopologyDihedral::atom3)
        .def_readwrite("atom4", &::pygcmc::model::TopologyDihedral::atom4)
        .def_readwrite("multiplicity", &::pygcmc::model::TopologyDihedral::multiplicity)
        .def_readwrite("angle", &::pygcmc::model::TopologyDihedral::angle)
        .def_readwrite("force_constant", &::pygcmc::model::TopologyDihedral::force_constant)
        .def_readwrite("improper", &::pygcmc::model::TopologyDihedral::improper)
        .def_readwrite("function_type", &::pygcmc::model::TopologyDihedral::function_type);

    // Bind TopologyDonor
    py::class_<::pygcmc::model::TopologyDonor>(model_module, "TopologyDonor")
        .def(py::init<>())
        .def_readwrite("donor_atom", &::pygcmc::model::TopologyDonor::donor_atom)
        .def_readwrite("hydrogen_atom", &::pygcmc::model::TopologyDonor::hydrogen_atom);

    // Bind TopologyAcceptor
    py::class_<::pygcmc::model::TopologyAcceptor>(model_module, "TopologyAcceptor")
        .def(py::init<>())
        .def_readwrite("acceptor_atom", &::pygcmc::model::TopologyAcceptor::acceptor_atom);

    // Bind TopologyCmap
    py::class_<::pygcmc::model::TopologyCmap>(model_module, "TopologyCmap")
        .def(py::init<>())
        .def_readwrite("atoms", &::pygcmc::model::TopologyCmap::atoms)
        .def_readwrite("function_type", &::pygcmc::model::TopologyCmap::function_type);

    // Bind StandardCmap
    py::class_<::pygcmc::model::StandardCmap>(model_module, "StandardCmap")
        .def(py::init<>())
        .def_readwrite("atoms", &::pygcmc::model::StandardCmap::atoms)
        .def_readwrite("raw_atoms", &::pygcmc::model::StandardCmap::raw_atoms)
        .def_readwrite("is_psf_format", &::pygcmc::model::StandardCmap::is_psf_format)
        .def_readwrite("function_type", &::pygcmc::model::StandardCmap::function_type);

    // Bind Topology class
    py::class_<::pygcmc::model::Topology>(model_module, "Topology")
        .def(py::init<>())
        .def("add_atom", &::pygcmc::model::Topology::add_atom)
        .def("add_bond", &::pygcmc::model::Topology::add_bond)
        .def("add_angle", &::pygcmc::model::Topology::add_angle)
        .def("add_dihedral", &::pygcmc::model::Topology::add_dihedral)
        .def("add_improper", &::pygcmc::model::Topology::add_improper)
        .def("add_donor", &::pygcmc::model::Topology::add_donor)
        .def("add_acceptor", &::pygcmc::model::Topology::add_acceptor)
        .def("add_nonbonded_exclusion", &::pygcmc::model::Topology::add_nonbonded_exclusion)
        .def("add_group", &::pygcmc::model::Topology::add_group)
        .def("add_cmap", static_cast<void (::pygcmc::model::Topology::*)(const std::array<int, 8>&)>(&::pygcmc::model::Topology::add_cmap), "Add CMAP in CHARMM format (8 atoms)")
        .def("add_cmap", static_cast<void (::pygcmc::model::Topology::*)(const std::array<int, 5>&, int)>(&::pygcmc::model::Topology::add_cmap), "Add CMAP in GROMACS format (5 atoms)", py::arg("atoms"), py::arg("function_type") = 1)
        .def("add_title", &::pygcmc::model::Topology::add_title)
        .def("get_num_atoms", &::pygcmc::model::Topology::get_num_atoms)
        .def("get_num_residues", &::pygcmc::model::Topology::get_num_residues)
        .def("get_num_segments", &::pygcmc::model::Topology::get_num_segments)
        .def("get_num_bonds", &::pygcmc::model::Topology::get_num_bonds)
        .def("get_num_angles", &::pygcmc::model::Topology::get_num_angles)
        .def("get_num_dihedrals", &::pygcmc::model::Topology::get_num_dihedrals)
        .def("get_num_impropers", &::pygcmc::model::Topology::get_num_impropers)
        .def("get_num_donors", &::pygcmc::model::Topology::get_num_donors)
        .def("get_num_acceptors", &::pygcmc::model::Topology::get_num_acceptors)
        .def("get_num_cmaps", &::pygcmc::model::Topology::get_num_cmaps)
        .def("get_num_groups", &::pygcmc::model::Topology::get_num_groups)
        .def("has_atom", &::pygcmc::model::Topology::has_atom)
        .def("has_residue", &::pygcmc::model::Topology::has_residue)
        .def("has_segment", &::pygcmc::model::Topology::has_segment)
        .def("has_bond", &::pygcmc::model::Topology::has_bond)
        .def("has_angle", &::pygcmc::model::Topology::has_angle)
        .def("has_dihedral", &::pygcmc::model::Topology::has_dihedral)
        .def("has_improper", &::pygcmc::model::Topology::has_improper)
        .def("has_donor", py::overload_cast<int, int>(&::pygcmc::model::Topology::has_donor, py::const_))
        .def("has_acceptor", &::pygcmc::model::Topology::has_acceptor)
        .def("has_cmap", py::overload_cast<>(&::pygcmc::model::Topology::has_cmap, py::const_))
        .def("has_cmap", py::overload_cast<const std::vector<int>&>(&::pygcmc::model::Topology::has_cmap, py::const_))
        .def("has_group", &::pygcmc::model::Topology::has_group)
        .def("get_atom", py::overload_cast<int>(&::pygcmc::model::Topology::get_atom, py::const_), py::return_value_policy::reference_internal, "Get atom by index")
        .def("get_residue", &::pygcmc::model::Topology::get_residue, py::return_value_policy::reference_internal, "Get residue by index")
        .def("get_segment", &::pygcmc::model::Topology::get_segment, py::return_value_policy::reference_internal, "Get segment by index")
        .def("get_group", &::pygcmc::model::Topology::get_group, py::return_value_policy::reference_internal, "Get group by index")
        .def("get_donor", &::pygcmc::model::Topology::get_donor, py::return_value_policy::reference_internal, "Get donor by index")
        .def("get_acceptor", &::pygcmc::model::Topology::get_acceptor, py::return_value_policy::reference_internal, "Get acceptor by index")
        .def("find_atom", &::pygcmc::model::Topology::find_atom)
        .def("find_residue", &::pygcmc::model::Topology::find_residue)
        .def("find_segment", &::pygcmc::model::Topology::find_segment)
        .def("get_cmaps", &::pygcmc::model::Topology::get_cmaps, py::return_value_policy::reference_internal, "Get all CMAP terms")
        .def_property_readonly("cmaps", [](const ::pygcmc::model::Topology& t) { return t.get_cmaps(); }, "Get all CMAP terms");

    // Bind TopologyResidue
    py::class_<::pygcmc::model::TopologyResidue>(model_module, "TopologyResidue")
        .def_readonly("name", &::pygcmc::model::TopologyResidue::name)
        .def_readonly("number", &::pygcmc::model::TopologyResidue::number)
        .def_readonly("atoms", &::pygcmc::model::TopologyResidue::atoms)
        .def_readonly("segment", &::pygcmc::model::TopologyResidue::segment);

    // Bind TopologyGroup
    py::class_<::pygcmc::model::TopologyGroup>(model_module, "TopologyGroup")
        .def_readonly("id", &::pygcmc::model::TopologyGroup::id)
        .def_readonly("atoms", &::pygcmc::model::TopologyGroup::atoms)
        .def_readonly("type", &::pygcmc::model::TopologyGroup::type);
}

} // namespace model
} // namespace bindings
} // namespace pygcmc