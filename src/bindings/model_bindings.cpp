#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "model/atom.hpp"
#include "model/residue.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {

void init_model(py::module& m) {
    // Create model submodule
    auto model = m.def_submodule("model", "Data model classes");
    
    // Bind Atom class
    py::class_<model::Atom, std::shared_ptr<model::Atom>>(model, "Atom")
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

    // Bind Residue class
    py::class_<model::Residue, std::shared_ptr<model::Residue>>(model, "Residue")
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
}

} // namespace bindings
} // namespace pygcmc 