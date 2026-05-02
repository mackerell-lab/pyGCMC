// src/bindings/platform/FragmentReservoirBindings.cpp
// Python bindings for FragmentReservoir

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../platform/cpu/movement/reservoir/fragment_reservoir.hpp"
#include "../../platform/cpu/movement/common/MovementUtils.hpp"

namespace py = pybind11;
using namespace pygcmc::platform::cpu::movement;

namespace pygcmc {
namespace bindings {
namespace platform {

void init_fragment_reservoir_bindings(py::module& m) {
    // Bind Quaternion class if not already bound
    if (!py::hasattr(m, "Quaternion")) {
        py::class_<Quaternion>(m, "Quaternion")
            .def(py::init<>())
            .def(py::init<double, double, double, double>())
            .def_readwrite("w", &Quaternion::w)
            .def_readwrite("x", &Quaternion::x)
            .def_readwrite("y", &Quaternion::y)
            .def_readwrite("z", &Quaternion::z)
            .def("normalize", &Quaternion::normalize);
    }

    // FragmentTemplate binding
    py::class_<FragmentTemplate>(m, "FragmentTemplate")
        .def(py::init<>())
        .def_readwrite("name", &FragmentTemplate::name)
        .def_readwrite("typeId", &FragmentTemplate::typeId)
        .def_readwrite("atoms", &FragmentTemplate::atoms)
        .def_readwrite("molecularWeight", &FragmentTemplate::molecularWeight)
        .def_readwrite("radius", &FragmentTemplate::radius)
        .def_readwrite("chemicalPotential", &FragmentTemplate::chemicalPotential)
        .def_readwrite("activity", &FragmentTemplate::activity)
        .def_readwrite("concentration", &FragmentTemplate::concentration)
        .def_readwrite("useCavityBias", &FragmentTemplate::useCavityBias)
        .def_readwrite("useConfigBias", &FragmentTemplate::useConfigBias)
        .def_readwrite("configBiasTrials", &FragmentTemplate::configBiasTrials)
        .def_readwrite("isRigid", &FragmentTemplate::isRigid)
        .def_readwrite("allowRotation", &FragmentTemplate::allowRotation)
        .def_readwrite("allowTranslation", &FragmentTemplate::allowTranslation)
        .def_readwrite("totalInsertions", &FragmentTemplate::totalInsertions)
        .def_readwrite("successfulInsertions", &FragmentTemplate::successfulInsertions)
        .def_readwrite("averageLifetime", &FragmentTemplate::averageLifetime)
        .def("getInsertionProbability", &FragmentTemplate::getInsertionProbability)
        .def("updateActivity", &FragmentTemplate::updateActivity)
        .def("calculateActivity", &FragmentTemplate::calculateActivity);

    // FragmentInstance binding
    py::class_<FragmentInstance>(m, "FragmentInstance")
        .def(py::init<>())
        .def_readwrite("instanceId", &FragmentInstance::instanceId)
        .def_readwrite("templateId", &FragmentInstance::templateId)
        .def_readwrite("residueIndex", &FragmentInstance::residueIndex)
        .def_readwrite("position", &FragmentInstance::position)
        .def_readwrite("centerOfMass", &FragmentInstance::centerOfMass)
        .def_readwrite("orientation", &FragmentInstance::orientation)
        .def_readwrite("velocity", &FragmentInstance::velocity)
        .def_readwrite("isActive", &FragmentInstance::isActive)
        .def_readwrite("isGhost", &FragmentInstance::isGhost)
        .def_readwrite("isFixed", &FragmentInstance::isFixed)
        .def_readwrite("energy_vdw", &FragmentInstance::energy_vdw)
        .def_readwrite("energy_elec", &FragmentInstance::energy_elec)
        .def_readwrite("energy_total", &FragmentInstance::energy_total)
        .def_readwrite("insertionTime", &FragmentInstance::insertionTime)
        .def_readwrite("lastMoveTime", &FragmentInstance::lastMoveTime)
        .def_readwrite("moveAttempts", &FragmentInstance::moveAttempts)
        .def_readwrite("acceptedMoves", &FragmentInstance::acceptedMoves)
        .def("getAcceptanceRate", &FragmentInstance::getAcceptanceRate)
        .def("getLifetime", &FragmentInstance::getLifetime)
        .def("needsNeighborUpdate", &FragmentInstance::needsNeighborUpdate);

    // FragmentReservoir::Config binding
    py::class_<FragmentReservoir::Config>(m, "FragmentReservoirConfig")
        .def(py::init<>())
        .def_readwrite("maxInstances", &FragmentReservoir::Config::maxInstances)
        .def_readwrite("maxGhosts", &FragmentReservoir::Config::maxGhosts)
        .def_readwrite("ghostRecycleRatio", &FragmentReservoir::Config::ghostRecycleRatio)
        .def_readwrite("autoCompact", &FragmentReservoir::Config::autoCompact)
        .def_readwrite("compactThreshold", &FragmentReservoir::Config::compactThreshold)
        .def_readwrite("trackStatistics", &FragmentReservoir::Config::trackStatistics)
        .def_readwrite("statisticsWindow", &FragmentReservoir::Config::statisticsWindow);

    // FragmentReservoir::Statistics binding
    py::class_<FragmentReservoir::Statistics>(m, "FragmentReservoirStatistics")
        .def(py::init<>())
        .def_readwrite("activeCountByType", &FragmentReservoir::Statistics::activeCountByType)
        .def_readwrite("ghostCountByType", &FragmentReservoir::Statistics::ghostCountByType)
        .def_readwrite("averageLifetimeByType", &FragmentReservoir::Statistics::averageLifetimeByType)
        .def_readwrite("acceptanceRateByType", &FragmentReservoir::Statistics::acceptanceRateByType)
        .def_readwrite("totalInsertions", &FragmentReservoir::Statistics::totalInsertions)
        .def_readwrite("totalDeletions", &FragmentReservoir::Statistics::totalDeletions)
        .def_readwrite("ghostRecycles", &FragmentReservoir::Statistics::ghostRecycles)
        .def_readwrite("memoryCompactions", &FragmentReservoir::Statistics::memoryCompactions)
        .def_readwrite("averageGhostLifetime", &FragmentReservoir::Statistics::averageGhostLifetime)
        .def_readwrite("peakActiveCount", &FragmentReservoir::Statistics::peakActiveCount)
        .def_readwrite("averageActiveCount", &FragmentReservoir::Statistics::averageActiveCount)
        .def_readwrite("insertionTimeMs", &FragmentReservoir::Statistics::insertionTimeMs)
        .def_readwrite("deletionTimeMs", &FragmentReservoir::Statistics::deletionTimeMs)
        .def_readwrite("queryTimeMs", &FragmentReservoir::Statistics::queryTimeMs)
        .def_readwrite("cacheHits", &FragmentReservoir::Statistics::cacheHits)
        .def_readwrite("cacheMisses", &FragmentReservoir::Statistics::cacheMisses)
        .def("print", &FragmentReservoir::Statistics::print)
        .def("reset", &FragmentReservoir::Statistics::reset);

    // FragmentReservoir binding
    py::class_<FragmentReservoir>(m, "FragmentReservoir")
        .def(py::init<>())
        .def(py::init<const FragmentReservoir::Config&>())

        // Template management
        .def("addTemplate", &FragmentReservoir::addTemplate,
             py::arg("template"),
             "Add a new fragment template")
        .def("loadTemplate", &FragmentReservoir::loadTemplate,
             py::arg("filename"),
             py::arg("name"),
             py::arg("chemicalPotential") = 0.0,
             "Load template from file")
        .def("getTemplate",
             py::overload_cast<int>(&FragmentReservoir::getTemplate),
             py::arg("templateId"),
             py::return_value_policy::reference,
             "Get template by ID")
        .def("getTemplate",
             py::overload_cast<const std::string&>(&FragmentReservoir::getTemplate),
             py::arg("name"),
             py::return_value_policy::reference,
             "Get template by name")
        .def("getTemplateCount", &FragmentReservoir::getTemplateCount,
             "Get number of templates")
        .def("updateChemicalPotential", &FragmentReservoir::updateChemicalPotential,
             py::arg("templateId"),
             py::arg("mu"),
             "Update chemical potential for a template")
        .def("updateActivity", &FragmentReservoir::updateActivity,
             py::arg("templateId"),
             py::arg("beta"),
             "Update activity for a template")
        .def("updateAllActivities", &FragmentReservoir::updateAllActivities,
             py::arg("beta"),
             "Update activities for all templates")

        // Instance management
        .def("createInstance", &FragmentReservoir::createInstance,
             py::arg("templateId"),
             py::arg("position"),
             py::arg("orientation") = Quaternion(),
             "Create a new fragment instance")
        .def("deleteInstance", &FragmentReservoir::deleteInstance,
             py::arg("instanceId"),
             "Delete an instance (convert to ghost)")
        .def("restoreInstance", &FragmentReservoir::restoreInstance,
             py::arg("instanceId"),
             py::arg("position"),
             py::arg("orientation"),
             "Restore a deleted (ghost) instance")
        .def("purgeInstance", &FragmentReservoir::purgeInstance,
             py::arg("instanceId"),
             "Permanently remove an instance")

        // Ghost management
        .def("recycleGhost", &FragmentReservoir::recycleGhost,
             py::arg("templateId"),
             "Recycle a ghost fragment")
        .def("purgeGhosts", &FragmentReservoir::purgeGhosts,
             py::arg("maxToKeep") = -1,
             "Purge old ghosts")
        .def("getGhostCount", &FragmentReservoir::getGhostCount,
             py::arg("templateId") = -1,
             "Get ghost count")

        // Query methods
        .def("getInstance",
             py::overload_cast<int>(&FragmentReservoir::getInstance),
             py::arg("instanceId"),
             py::return_value_policy::reference_internal,
             "Get instance by ID")
        .def("getActiveInstances", &FragmentReservoir::getActiveInstances,
             py::arg("templateId") = -1,
             "Get active instances")
        .def("getActiveCount", &FragmentReservoir::getActiveCount,
             py::arg("templateId") = -1,
             "Get active instance count")

        // Statistics
        .def("getStatistics", &FragmentReservoir::getStatistics,
             py::return_value_policy::copy,
             "Get statistics")
        .def("resetStatistics", &FragmentReservoir::resetStatistics,
             "Reset statistics")
        .def("printStatistics", &FragmentReservoir::printStatistics,
             "Print statistics")

        // Memory management
        .def("compact", &FragmentReservoir::compact,
             "Compact memory")
        .def("getFragmentation", &FragmentReservoir::getFragmentation,
             "Get memory fragmentation ratio");
}

} // namespace platform
} // namespace bindings
} // namespace pygcmc
