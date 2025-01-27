#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "model/atom.hpp"
#include "model/residue.hpp"
#include "model/topology.hpp"

namespace py = pybind11;

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
        .def("centerOfMass", &model::Residue::centerOfMass)
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
        .def("add_cmap", &model::Topology::add_cmap)
        .def("add_title", &model::Topology::add_title)
        .def("get_num_atoms", &model::Topology::get_num_atoms)
        .def("get_num_residues", &model::Topology::get_num_residues)
        .def("get_num_segments", &model::Topology::get_num_segments)
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
}

} // namespace bindings
} // namespace pygcmc 