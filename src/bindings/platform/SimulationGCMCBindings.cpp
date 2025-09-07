// GCMCBindings.cpp - Python bindings for GCMC simulation

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
// Use new EnergyAPI instead of simulation.hpp
#include "../../platform/cpu/energy/EnergyAPI.hpp"
#include "../../platform/cpu/energy/common/EnergyInterface.hpp"
// GCMC specific headers
#include "../../platform/cpu/movement/gcmc/GCMCEngine.hpp"
#include "../../platform/cpu/movement/gcmc/GCMCAcceptance.hpp"
#include "../../platform/cpu/movement/reservoir/fragment_reservoir.hpp"
#include "../../platform/cpu/movement/bias/CavityBias.hpp"
#include "../../model/montecarlo/MCStructures.hpp"

namespace py = pybind11;

// Type aliases for clarity
using MCState = ::pygcmc::model::montecarlo::MCState;
using FragmentTemplate = ::pygcmc::platform::cpu::movement::FragmentTemplate;
using FragmentReservoir = ::pygcmc::platform::cpu::movement::FragmentReservoir;
using FragmentInstance = ::pygcmc::platform::cpu::movement::FragmentInstance;
using GCMCEngine = ::pygcmc::platform::cpu::movement::gcmc::GCMCEngine;
using GCMCAcceptance = ::pygcmc::platform::cpu::movement::gcmc::GCMCAcceptance;
using CavityManager = ::pygcmc::platform::cpu::movement::CavityManager;
using Vector3 = ::pygcmc::platform::cpu::movement::Vector3;
using Quaternion = ::pygcmc::platform::cpu::movement::Quaternion;

namespace pygcmc {
namespace bindings {
namespace platform {

void init_gcmc_bindings(py::module& m) {
    // Energy method enum
    py::enum_<::pygcmc::platform::cpu::EnergyMethod>(m, "EnergyMethod")
        .value("DIRECT", ::pygcmc::platform::cpu::EnergyMethod::DIRECT)
        .value("EWALD", ::pygcmc::platform::cpu::EnergyMethod::EWALD)
        .value("PME", ::pygcmc::platform::cpu::EnergyMethod::PME);
    
    // Note: FragmentTemplate, FragmentInstance, and FragmentReservoir are already
    // registered in FragmentReservoirBindings.cpp, so we don't re-register them here.
    
    // GCMCEngine class
    py::class_<GCMCEngine>(m, "GCMCEngine")
        .def(py::init<>())
        .def("initialize", &GCMCEngine::initialize,
             py::arg("state"), py::arg("reservoir"),
             "Initialize the GCMC engine with state and reservoir")
        .def("setSeed", &GCMCEngine::setSeed,
             py::arg("seed"), "Set random seed")
        .def("setTemperature", &GCMCEngine::setTemperature,
             py::arg("temperature"), "Set temperature in Kelvin")
        .def("setCutoff", &GCMCEngine::setCutoff,
             py::arg("cutoff"), "Set cutoff distance in Angstroms")
        .def("attemptInsertion", &GCMCEngine::attemptInsertion,
             py::arg("typeId"), "Attempt an insertion move")
        .def("attemptDeletion", &GCMCEngine::attemptDeletion,
             py::arg("typeId"), "Attempt a deletion move")
        .def("attemptTranslation", &GCMCEngine::attemptTranslation,
             py::arg("residueIdx"), "Attempt a translation move")
        .def("attemptRotation", &GCMCEngine::attemptRotation,
             py::arg("residueIdx"), "Attempt a rotation move")
        .def("setAcceptanceCalculator", &GCMCEngine::setAcceptanceCalculator,
             py::arg("calculator"), py::keep_alive<1, 2>(),
             "Set the acceptance calculator")
        .def("setCavityManager", &GCMCEngine::setCavityManager,
             py::arg("manager"), py::keep_alive<1, 2>(),
             "Set the cavity manager")
        .def("getAcceptanceRate", &GCMCEngine::getAcceptanceRate,
             "Get the overall acceptance rate")
        .def("synchronizeStateWithReservoir", &GCMCEngine::synchronizeStateWithReservoir,
             "Synchronize MCState with the reservoir's active fragments");
    
    // MoveResult struct for GCMCEngine
    py::class_<GCMCEngine::MoveResult>(m, "GCMCMoveResult")
        .def_readonly("accepted", &GCMCEngine::MoveResult::accepted)
        .def_readonly("deltaE", &GCMCEngine::MoveResult::deltaE)
        .def_readonly("energyBefore", &GCMCEngine::MoveResult::energyBefore)
        .def_readonly("energyAfter", &GCMCEngine::MoveResult::energyAfter)
        .def_readonly("bias", &GCMCEngine::MoveResult::bias)
        .def_readonly("acceptanceProbability", &GCMCEngine::MoveResult::acceptanceProbability)
        .def_readonly("residueIndex", &GCMCEngine::MoveResult::residueIndex);
    
    // GCMCAcceptance class
    py::class_<GCMCAcceptance>(m, "GCMCAcceptance")
        .def(py::init<>())
        .def("setTemperature", &GCMCAcceptance::setTemperature,
             py::arg("temperature"), "Set temperature in Kelvin")
        .def("setVolume", &GCMCAcceptance::setVolume,
             py::arg("volume"), "Set volume in nm^3")
        .def("setPressure", &GCMCAcceptance::setPressure,
             py::arg("pressure"), "Set pressure in bar")
        .def("setChemicalPotential", &GCMCAcceptance::setChemicalPotential,
             py::arg("typeId"), py::arg("mu"),
             "Set chemical potential for a type")
        .def("setActivity", &GCMCAcceptance::setActivity,
             py::arg("typeId"), py::arg("activity"),
             "Set activity for a type")
        .def("calculateInsertionProbability", 
             &GCMCAcceptance::calculateInsertionProbability,
             py::arg("typeId"), py::arg("currentNumber"), 
             py::arg("deltaE"), py::arg("bias"),
             "Calculate insertion acceptance probability")
        .def("calculateDeletionProbability",
             &GCMCAcceptance::calculateDeletionProbability,
             py::arg("typeId"), py::arg("currentNumber"),
             py::arg("deltaE"), py::arg("bias"),
             "Calculate deletion acceptance probability")
        .def("calculateTranslationProbability",
             &GCMCAcceptance::calculateTranslationProbability,
             py::arg("deltaE"), py::arg("bias") = 1.0,
             "Calculate translation acceptance probability")
        .def("setSeed", &GCMCAcceptance::setSeed,
             py::arg("seed"), "Set random seed for acceptance decisions")
        .def("acceptMove", &GCMCAcceptance::acceptMove,
             py::arg("probability"), "Accept or reject based on probability");
    
    // CavityManager class
    py::class_<CavityManager>(m, "CavityManager")
        .def(py::init<double, double>(),
             py::arg("gridSpacing") = 2.0, py::arg("probeRadius") = 1.4,
             "Create cavity manager with grid spacing and probe radius")
        .def("findCavities", &CavityManager::findCavities,
             py::arg("state"), "Find cavities in the system")
        .def("invalidateCache", &CavityManager::invalidateCache,
             "Invalidate the cavity cache")
        .def("getCavityVolume", &CavityManager::getCavityVolume,
             py::arg("state"), "Get total cavity volume in nm^3")
        .def("getCavityVolumeFraction", &CavityManager::getCavityVolumeFraction,
             py::arg("state"), "Get cavity volume fraction")
        .def("setGridSpacing", &CavityManager::setGridSpacing,
             py::arg("spacing"), "Set grid spacing")
        .def("setProbeRadius", &CavityManager::setProbeRadius,
             py::arg("radius"), "Set probe radius")
        .def("getCavityCount", &CavityManager::getCavityCount,
             "Get number of cavities found");
}

} // namespace platform
} // namespace bindings
} // namespace pygcmc