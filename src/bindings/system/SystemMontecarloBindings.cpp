// src/bindings/system/SystemMontecarloBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "system/SystemModule.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace system {

void init_montecarlo_bindings(py::module& m, py::module&) {
    // Bind MonteCarloSystem
    py::class_<pygcmc::system::MonteCarloSystem>(m, "MonteCarloSystem")
        .def(py::init<>())
        .def("initialize", &pygcmc::system::MonteCarloSystem::initialize)
        .def("set_force_field", &pygcmc::system::MonteCarloSystem::setForceField)
        .def("initialize_force_field", &pygcmc::system::MonteCarloSystem::initializeForceField)
        .def("initialize_from_molecular", [](pygcmc::system::MonteCarloSystem& self, py::object molecular) {
            if (molecular.is_none()) {
                throw py::value_error("Molecular object cannot be None");
            }

            try {
                // First try MolecularSystem
                auto* molSys = molecular.cast<pygcmc::system::MolecularSystem*>();
                if (molSys) {
                    self.initializeFromMolecular(molSys->get_molecular());
                    return;
                }
            } catch (py::cast_error&) {}

            try {
                // Then try Molecular directly
                auto mol = molecular.cast<std::shared_ptr<pygcmc::model::Molecular>>();
                if (mol) {
                    self.initializeFromMolecular(mol);
                    return;
                }
            } catch (py::cast_error&) {}

            throw py::type_error("Argument must be either MolecularSystem or Molecular");
        })
        .def("add_movement_molecules", [](pygcmc::system::MonteCarloSystem& self, py::list molecules) {
            std::vector<pygcmc::system::montecarlo::MovementMolecularInfo> mol_vec;
            for (const auto& mol : molecules) {
                mol_vec.push_back(mol.cast<pygcmc::system::montecarlo::MovementMolecularInfo>());
            }
            self.addMovementMolecules(mol_vec);
        }, py::arg("molecules"), "Add movement molecules for GCMC simulation")
        .def("get_type_maps", &pygcmc::system::MonteCarloSystem::getTypeMaps, py::return_value_policy::reference)
        .def("insert_residue", &pygcmc::system::MonteCarloSystem::insertResidue)
        .def("remove_residue", &pygcmc::system::MonteCarloSystem::removeResidue)
        .def("translate_residue", &pygcmc::system::MonteCarloSystem::translateResidue)
        .def("calc_non_bonded_energy", &pygcmc::system::MonteCarloSystem::calcNonBondedEnergy)
        .def("calc_total_energy", &pygcmc::system::MonteCarloSystem::calcTotalEnergy)
        .def("get_state", (const pygcmc::model::MCState& (pygcmc::system::MonteCarloSystem::*)() const) &pygcmc::system::MonteCarloSystem::getState, py::return_value_policy::reference)
        .def("get_state_mutable", (pygcmc::model::MCState& (pygcmc::system::MonteCarloSystem::*)()) &pygcmc::system::MonteCarloSystem::getState, py::return_value_policy::reference)
        .def("get_active_atom_count", &pygcmc::system::MonteCarloSystem::getActiveAtomCount)
        .def("get_active_residue_count", &pygcmc::system::MonteCarloSystem::getActiveResidueCount)
        .def("set_switching_function", &pygcmc::system::MonteCarloSystem::setSwitchingFunction,
             py::arg("enable"), py::arg("r_on") = 1.0f, py::arg("r_off") = 1.2f,
             "Set or disable CHARMM-style smooth switching function")
        .def("calculate_switching_function", &pygcmc::system::MonteCarloSystem::calculateSwitchingFunction,
             py::arg("r"), "Calculate switching function value at the given distance")
        .def("is_using_switching_function", &pygcmc::system::MonteCarloSystem::isUsingSwitchingFunction,
             "Get whether the switching function is currently enabled")
        .def("get_switching_r_on", &pygcmc::system::MonteCarloSystem::getSwitchingROn,
             "Get the inner cutoff radius")
        .def("get_switching_r_off", &pygcmc::system::MonteCarloSystem::getSwitchingROff,
             "Get the outer cutoff radius")
        .def("apply_switching_to_state", &pygcmc::system::MonteCarloSystem::applySwitchingToState,
             py::arg("state"), "Apply the current switching function settings to an external state object");

    // Bind MovementMolecularInfo
    py::class_<pygcmc::system::montecarlo::MovementMolecularInfo>(m, "MovementMolecularInfo")
        .def(py::init<std::shared_ptr<pygcmc::model::Molecular>, int>(),
             py::arg("molecular"),
             py::arg("maxCopies"))
        .def_readwrite("molecular", &pygcmc::system::montecarlo::MovementMolecularInfo::molecular)
        .def_readwrite("maxCopies", &pygcmc::system::montecarlo::MovementMolecularInfo::maxCopies);
}

} // namespace system
} // namespace bindings
} // namespace pygcmc
