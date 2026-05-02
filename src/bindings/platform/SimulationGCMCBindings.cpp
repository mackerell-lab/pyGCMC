// GCMCBindings.cpp - Python bindings for GCMC simulation

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <utility>
// Use new EnergyAPI instead of simulation.hpp
#include "../../platform/cpu/energy/EnergyAPI.hpp"
#include "../../platform/cpu/energy/common/EnergyInterface.hpp"
// GCMC specific headers
#include "../../platform/cpu/movement/gcmc/GCMCEngine.hpp"
#include "../../platform/cpu/movement/gcmc/GCMCAcceptance.hpp"
#include "../../platform/cpu/movement/reservoir/fragment_reservoir.hpp"
#include "../../platform/cpu/movement/bias/CavityBias.hpp"
#include "../../platform/cpu/simulation/impl/GCMCSimulation.hpp"
#include "../../model/montecarlo/MCStructures.hpp"

namespace py = pybind11;

// Type aliases for clarity
using MCState = ::pygcmc::model::montecarlo::MCState;
using FragmentTemplate = ::pygcmc::platform::cpu::movement::FragmentTemplate;
using FragmentReservoir = ::pygcmc::platform::cpu::movement::FragmentReservoir;
using FragmentInstance = ::pygcmc::platform::cpu::movement::FragmentInstance;
using GCMCEngine = ::pygcmc::platform::cpu::movement::gcmc::GCMCEngine;
using GCMCAcceptance = ::pygcmc::platform::cpu::movement::gcmc::GCMCAcceptance;
using GCMCStatistics = ::pygcmc::platform::cpu::movement::gcmc::GCMCStatistics;
using CavityManager = ::pygcmc::platform::cpu::movement::CavityManager;
using Vector3 = ::pygcmc::platform::cpu::movement::Vector3;
using Quaternion = ::pygcmc::platform::cpu::movement::Quaternion;
using GCMCCPUSimulation = ::pygcmc::platform::cpu::simulation::GCMCSimulation;

namespace {

py::dict cpu_stats_to_dict(const GCMCCPUSimulation::Statistics& stats) {
    py::dict result;
    result["totalSteps"] = stats.totalSteps;
    result["acceptedMoves"] = stats.acceptedMoves;
    result["acceptanceRate"] = stats.acceptanceRate;
    result["moveAttempts"] = stats.moveAttempts;
    result["moveAccepted"] = stats.moveAccepted;
    result["moveAcceptanceRates"] = stats.moveAcceptanceRates;
    result["fragmentCounts"] = stats.fragmentCounts;
    result["fragmentDensities"] = stats.fragmentDensities;
    result["fragmentAcceptanceRates"] = stats.fragmentAcceptanceRates;
    result["currentEnergy"] = stats.currentEnergy;
    result["averageEnergy"] = stats.averageEnergy;
    result["energyStdDev"] = stats.energyStdDev;
    result["energyHistory"] = stats.energyHistory;
    result["totalTime"] = stats.totalTime;
    result["timePerStep"] = stats.timePerStep;
    result["stepsPerSecond"] = stats.stepsPerSecond;
    return result;
}

py::dict cpu_counters_to_dict(const GCMCCPUSimulation::Counters& counters) {
    py::dict result;
    result["attemptsInsert"] = counters.attemptsInsert;
    result["acceptsInsert"] = counters.acceptsInsert;
    result["attemptsDelete"] = counters.attemptsDelete;
    result["acceptsDelete"] = counters.acceptsDelete;
    result["attemptsTranslate"] = counters.attemptsTranslate;
    result["acceptsTranslate"] = counters.acceptsTranslate;
    result["attemptsRotate"] = counters.attemptsRotate;
    result["acceptsRotate"] = counters.acceptsRotate;
    result["attemptsSinceLastPrint"] = counters.attemptsSinceLastPrint;
    result["acceptsSinceLastPrint"] = counters.acceptsSinceLastPrint;
    result["insDelOverallRate"] = counters.getInsDelOverallRate();
    return result;
}

py::dict cpu_fragment_to_dict(const GCMCCPUSimulation::FragmentInfo& fragment) {
    py::dict result;
    result["name"] = fragment.name;
    result["typeId"] = fragment.typeId;
    result["concentration"] = fragment.concentration;
    result["chemicalPotential"] = fragment.chemicalPotential;
    result["activity"] = fragment.activity;
    result["probability"] = fragment.probability;
    result["maxCount"] = fragment.maxCount;
    result["currentCount"] = fragment.currentCount;
    result["confBiasTrials"] = fragment.confBiasTrials;
    result["insertAttempts"] = fragment.insertAttempts;
    result["insertAccepted"] = fragment.insertAccepted;
    result["deleteAttempts"] = fragment.deleteAttempts;
    result["deleteAccepted"] = fragment.deleteAccepted;
    return result;
}

py::list cpu_fragments_to_list(const std::vector<GCMCCPUSimulation::FragmentInfo>& fragments) {
    py::list result;
    for (const auto& fragment : fragments) {
        result.append(cpu_fragment_to_dict(fragment));
    }
    return result;
}

std::string cpu_move_type_to_string(GCMCCPUSimulation::AcceptanceRecord::MoveType moveType) {
    switch (moveType) {
        case GCMCCPUSimulation::AcceptanceRecord::INSERT:
            return "insertion";
        case GCMCCPUSimulation::AcceptanceRecord::DELETE:
            return "deletion";
        case GCMCCPUSimulation::AcceptanceRecord::TRANSLATE:
            return "translation";
        case GCMCCPUSimulation::AcceptanceRecord::ROTATE:
            return "rotation";
    }
    return "unknown";
}

py::dict cpu_acceptance_record_to_dict(const GCMCCPUSimulation::AcceptanceRecord& record) {
    py::dict result;
    result["move"] = cpu_move_type_to_string(record.moveType);
    result["requestedMove"] = cpu_move_type_to_string(record.requestedMoveType);
    result["species"] = record.species;
    result["nBefore"] = record.nBefore;
    result["cbmcTrials"] = record.cbmcTrials;
    result["step"] = record.step;
    result["beta"] = record.beta;
    result["deltaU"] = record.deltaU;
    result["betaDeltaU"] = record.betaDeltaU;
    result["mu"] = record.mu;
    result["betaMu"] = record.betaMu;
    result["z"] = record.z;
    result["qForward"] = record.qForward;
    result["qReverse"] = record.qReverse;
    result["proposalRatio"] = record.proposalRatio;
    result["vEff"] = record.vEff;
    result["vBox"] = record.vBox;
    result["cavityFraction"] = record.cavityFraction;
    result["rosenbluthWeight"] = record.rosenbluthWeight;
    result["cbmcSelectedEnergy"] = record.cbmcSelectedEnergy;
    result["cbmcLogWOverK"] = record.cbmcLogWOverK;
    result["cbmcTrialEnergies"] = record.cbmcTrialEnergies;
    result["pAcc"] = record.pAcc;
    result["bias"] = record.bias;
    result["u"] = record.u;
    result["accepted"] = record.accepted;
    result["wForward"] = record.wForward;
    result["wReverse"] = record.wReverse;
    result["wCavity"] = record.wCavity;
    return result;
}

py::list cpu_acceptance_records_to_list(const std::vector<GCMCCPUSimulation::AcceptanceRecord>& records) {
    py::list result;
    for (const auto& record : records) {
        result.append(cpu_acceptance_record_to_dict(record));
    }
    return result;
}

py::dict run_gcmc_cpu_config(
    const GCMCCPUSimulation::Config& config,
    const std::string& resumeCheckpoint,
    const std::string& dumpAccept,
    const std::string& dumpParams,
    size_t diagnosticsBuffer) {
    GCMCCPUSimulation sim(config);

    bool initialized = false;
    {
        py::gil_scoped_release release;
        initialized = sim.initialize();
        if (!initialized && !dumpParams.empty()) {
            sim.dumpParamsJson(dumpParams);
        }
    }

    py::dict result;
    result["initialized"] = initialized;
    if (!initialized) {
        result["returncode"] = 2;
        result["ran"] = false;
        result["finalized"] = false;
        return result;
    }

    if (!dumpAccept.empty()) {
        sim.enableDiagnostics(diagnosticsBuffer);
    }

    bool checkpointLoaded = true;
    bool ran = false;
    {
        py::gil_scoped_release release;
        if (!dumpParams.empty()) {
            sim.dumpParamsJson(dumpParams);
        }
        if (!resumeCheckpoint.empty()) {
            checkpointLoaded = sim.loadCheckpoint(resumeCheckpoint);
        }
        if (checkpointLoaded) {
            ran = sim.run();
            if (ran) {
                sim.finalize();
                if (!dumpAccept.empty()) {
                    sim.dumpAcceptanceLog(dumpAccept);
                }
            }
        }
    }

    result["checkpointLoaded"] = checkpointLoaded;
    result["ran"] = ran;
    result["finalized"] = ran;
    result["returncode"] = checkpointLoaded ? (ran ? 0 : 3) : 4;
    result["statistics"] = cpu_stats_to_dict(sim.getStatistics());
    result["counters"] = cpu_counters_to_dict(sim.getCounters());
    result["fragments"] = cpu_fragments_to_list(sim.getFragmentInfo());
    return result;
}

} // namespace

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
             py::keep_alive<1, 2>(),  // Keep MCState alive as long as the engine exists
             py::keep_alive<1, 3>(),  // Keep FragmentReservoir alive as long as the engine exists
             "Initialize the GCMC engine with state and reservoir")
        .def("setSeed", &GCMCEngine::setSeed,
             py::arg("seed"), "Set random seed")
        .def("setTemperature", &GCMCEngine::setTemperature,
             py::arg("temperature"), "Set temperature in Kelvin")
        .def("setCutoff", &GCMCEngine::setCutoff,
             py::arg("cutoff"), "Set cutoff distance in Angstroms")
        .def(
            "setRegionConstraintFromSpec",
            [](GCMCEngine& self,
               const std::string& regionSpec,
               double boxX,
               double boxY,
               double boxZ) {
                using RegionConstraint = ::pygcmc::platform::cpu::movement::RegionConstraint;
                auto constraint = RegionConstraint::parseRegion(regionSpec, Vector3(boxX, boxY, boxZ));
                self.setRegionConstraint(std::move(constraint));
            },
            py::arg("regionSpec"),
            py::arg("boxX"),
            py::arg("boxY"),
            py::arg("boxZ"),
            "Set a region constraint from a gcmc_region spec (nm units)")
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
        .def("setCBMCTrialsPerType", &GCMCEngine::setCBMCTrialsPerType,
             py::arg("trials"),
             "Set CBMC trial numbers per fragment type")
        .def("getAcceptanceRate", &GCMCEngine::getAcceptanceRate,
             "Get the overall acceptance rate")
        .def("synchronizeStateWithReservoir", &GCMCEngine::synchronizeStateWithReservoir,
             "Synchronize MCState with the reservoir's active fragments")
        // Dynamic configuration
        .def("setConfigValue", &GCMCEngine::setConfigValue,
             py::arg("key"), py::arg("value"),
             "Set a configuration value dynamically")
        .def("getConfigValue", &GCMCEngine::getConfigValue,
             py::arg("key"),
             "Get a configuration value")
        // Statistics management
        .def("enableStatistics", &GCMCEngine::enableStatistics,
             py::arg("enable"),
             "Enable or disable statistics collection")
        .def("setStatisticsInterval", &GCMCEngine::setStatisticsInterval,
             py::arg("interval"),
             "Set the interval for statistics sampling")
        .def("getStatistics",
             (GCMCStatistics& (GCMCEngine::*)()) &GCMCEngine::getStatistics,
             py::return_value_policy::reference_internal,
             "Get the statistics collector")
        // New methods for residue position and orientation queries
        .def("getResiduePosition", &GCMCEngine::getResiduePosition,
             py::arg("residueIdx"), "Get position of residue")
        .def("getResidueOrientation", &GCMCEngine::getResidueOrientation,
             py::arg("residueIdx"), "Get orientation of residue");

    // GCMCStatistics class
    py::class_<GCMCStatistics>(m, "GCMCStatistics")
        .def(py::init<>())
        .def("setAutoAdjust", &GCMCStatistics::setAutoAdjust,
             py::arg("enable"), "Enable/disable auto-adjustment of sampling frequency")
        .def("setSamplingInterval", &GCMCStatistics::setSamplingInterval,
             py::arg("interval"), "Set sampling interval")
        .def("getSamplingInterval", &GCMCStatistics::getSamplingInterval,
             "Get current sampling interval")
        .def("shouldSample", &GCMCStatistics::shouldSample,
             py::arg("step"), "Check if should sample at this step")
        .def("addSample", &GCMCStatistics::addSample,
             py::arg("step"), py::arg("particleCount"), py::arg("energy"),
             py::arg("acceptanceRate"), py::arg("temperature"), py::arg("activity"),
             "Add a sample")
        .def("getParticleStats", &GCMCStatistics::getParticleStats,
             "Get particle count statistics")
        .def("getEnergyStats", &GCMCStatistics::getEnergyStats,
             "Get energy statistics")
        .def("getCurrentVariance", &GCMCStatistics::getCurrentVariance,
             "Get current variance for auto-adjustment")
        .def("clear", &GCMCStatistics::clear,
             "Clear all samples");

    // Stats structure
    py::class_<GCMCStatistics::Stats>(m, "GCMCStats")
        .def_readonly("mean", &GCMCStatistics::Stats::mean)
        .def_readonly("variance", &GCMCStatistics::Stats::variance)
        .def_readonly("min", &GCMCStatistics::Stats::min)
        .def_readonly("max", &GCMCStatistics::Stats::max)
        .def_readonly("count", &GCMCStatistics::Stats::count);

    // MoveResult struct for GCMCEngine
    py::class_<GCMCEngine::MoveResult>(m, "GCMCMoveResult")
        .def_readonly("accepted", &GCMCEngine::MoveResult::accepted)
        .def_readonly("deltaE", &GCMCEngine::MoveResult::deltaE)
        .def_readonly("energyBefore", &GCMCEngine::MoveResult::energyBefore)
        .def_readonly("energyAfter", &GCMCEngine::MoveResult::energyAfter)
        .def_readonly("bias", &GCMCEngine::MoveResult::bias)
        .def_readonly("acceptanceProbability", &GCMCEngine::MoveResult::acceptanceProbability)
        .def_readonly("residueIndex", &GCMCEngine::MoveResult::residueIndex)
        .def_readonly("position", &GCMCEngine::MoveResult::position)
        .def_readonly("rosenbluthWeight", &GCMCEngine::MoveResult::rosenbluthWeight)
        .def_readonly("cavityBiasComponent", &GCMCEngine::MoveResult::cavityBiasComponent)
        .def_readonly("effectiveVolume", &GCMCEngine::MoveResult::effectiveVolume)
        .def_readonly("cbmcTrialsUsed", &GCMCEngine::MoveResult::cbmcTrialsUsed);

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
        .def("setThermalLambda", &GCMCAcceptance::setThermalLambda,
             py::arg("typeId"), py::arg("lambdaNm"),
             "Set thermal de Broglie wavelength (nm) for a type")
        .def("getThermalLambda", &GCMCAcceptance::getThermalLambda,
             py::arg("typeId"),
             "Get thermal de Broglie wavelength (nm) for a type")
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
        .def(
            "calculate_insertion_probability_detailed",
            [](GCMCAcceptance& self,
               int typeId,
               int currentNumber,
               double deltaE,
               double cavityFraction,
               double lambdaNm,
               double rosenbluthWeight,
               int cbmcTrials,
               double proposalLogRatio) {
                GCMCAcceptance::GrandCanonicalInsertionTerms terms;
                terms.typeId = typeId;
                terms.countBefore = currentNumber;
                terms.deltaE = deltaE;
               terms.cavityFraction = cavityFraction;
               terms.lambdaNm = lambdaNm;
               terms.rosenbluthWeight = rosenbluthWeight;
               terms.cbmcTrials = cbmcTrials;
               terms.proposalLogRatio = proposalLogRatio;
                double logRatio = 0.0;
                double probability = self.calculateInsertionProbabilityDetailed(terms, &logRatio);
                py::dict result;
                result["probability"] = probability;
                result["logRatio"] = logRatio;
                return result;
            },
            py::arg("typeId"),
            py::arg("currentNumber"),
            py::arg("deltaE"),
            py::arg("cavityFraction"),
            py::arg("lambdaNm") = 1.0,
            py::arg("rosenbluthWeight") = 1.0,
            py::arg("cbmcTrials") = 1,
            py::arg("proposalLogRatio") = 0.0,
            "Calculate insertion probability with detailed balance terms")
        .def(
            "calculate_deletion_probability_detailed",
            [](GCMCAcceptance& self,
               int typeId,
               int currentNumber,
               double deltaE,
               double cavityFraction,
               double lambdaNm,
               double rosenbluthWeight,
               int cbmcTrials,
               double proposalLogRatio) {
                GCMCAcceptance::GrandCanonicalDeletionTerms terms;
                terms.typeId = typeId;
                terms.countBefore = currentNumber;
                terms.deltaE = deltaE;
               terms.cavityFraction = cavityFraction;
               terms.lambdaNm = lambdaNm;
               terms.rosenbluthWeight = rosenbluthWeight;
               terms.cbmcTrials = cbmcTrials;
               terms.proposalLogRatio = proposalLogRatio;
                double logRatio = 0.0;
                double probability = self.calculateDeletionProbabilityDetailed(terms, &logRatio);
                py::dict result;
                result["probability"] = probability;
                result["logRatio"] = logRatio;
                return result;
            },
            py::arg("typeId"),
            py::arg("currentNumber"),
            py::arg("deltaE"),
            py::arg("cavityFraction"),
            py::arg("lambdaNm") = 1.0,
            py::arg("rosenbluthWeight") = 1.0,
            py::arg("cbmcTrials") = 1,
            py::arg("proposalLogRatio") = 0.0,
            "Calculate deletion probability with detailed balance terms")
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
             py::arg("state"), py::arg("speciesId") = -1,
             "Find cavities in the system (optionally for a species)")
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

    py::class_<GCMCCPUSimulation::Config>(m, "GCMCCPUConfig")
        .def(py::init<>())
        .def_readwrite("inputFile", &GCMCCPUSimulation::Config::inputFile)
        .def_readwrite("outputPrefix", &GCMCCPUSimulation::Config::outputPrefix)
        .def_readwrite("printFrequency", &GCMCCPUSimulation::Config::printFrequency)
        .def_readwrite("trajectoryFrequency", &GCMCCPUSimulation::Config::trajectoryFrequency)
        .def_readwrite("checkpointFrequency", &GCMCCPUSimulation::Config::checkpointFrequency)
        .def_readwrite("verbose", &GCMCCPUSimulation::Config::verbose)
        .def_readwrite("randomSeed", &GCMCCPUSimulation::Config::randomSeed)
        .def_readwrite("strictInpKeys", &GCMCCPUSimulation::Config::strictInpKeys)
        .def_readwrite("strictInpWarnings", &GCMCCPUSimulation::Config::strictInpWarnings)
        .def_readwrite("enableStatistics", &GCMCCPUSimulation::Config::enableStatistics)
        .def_readwrite("movesPerStep", &GCMCCPUSimulation::Config::movesPerStep)
        .def_readwrite("statisticsInterval", &GCMCCPUSimulation::Config::statisticsInterval)
        .def_readwrite("storeProbabilities", &GCMCCPUSimulation::Config::storeProbabilities)
        .def_readwrite("enableAdaptiveSampling", &GCMCCPUSimulation::Config::enableAdaptiveSampling)
        .def_readwrite("enableEnergyMinimization", &GCMCCPUSimulation::Config::enableEnergyMinimization)
        .def_readwrite("convergenceTolerance", &GCMCCPUSimulation::Config::convergenceTolerance)
        .def_readwrite("maxMoleculesPerType", &GCMCCPUSimulation::Config::maxMoleculesPerType);

    py::class_<GCMCCPUSimulation>(m, "GCMCCPUSimulation")
        .def(py::init<const GCMCCPUSimulation::Config&>(), py::arg("config"))
        .def("initialize", &GCMCCPUSimulation::initialize,
             py::call_guard<py::gil_scoped_release>())
        .def("run", &GCMCCPUSimulation::run,
             py::call_guard<py::gil_scoped_release>())
        .def("finalize", &GCMCCPUSimulation::finalize,
             py::call_guard<py::gil_scoped_release>())
        .def("stop", &GCMCCPUSimulation::stop)
        .def("is_running", &GCMCCPUSimulation::isRunning)
        .def("save_trajectory", &GCMCCPUSimulation::saveTrajectory,
             py::arg("filename"),
             py::call_guard<py::gil_scoped_release>())
        .def("save_topology", &GCMCCPUSimulation::saveTopology,
             py::arg("filename"),
             py::call_guard<py::gil_scoped_release>())
        .def("save_checkpoint", &GCMCCPUSimulation::saveCheckpoint,
             py::arg("filename"),
             py::call_guard<py::gil_scoped_release>())
        .def("load_checkpoint", &GCMCCPUSimulation::loadCheckpoint,
             py::arg("filename"),
             py::call_guard<py::gil_scoped_release>())
        .def("enable_diagnostics", &GCMCCPUSimulation::enableDiagnostics,
             py::arg("bufferSize") = 4096)
        .def("is_diagnostics_enabled", &GCMCCPUSimulation::isDiagnosticsEnabled)
        .def("dump_acceptance_log", &GCMCCPUSimulation::dumpAcceptanceLog,
             py::arg("filename"),
             py::call_guard<py::gil_scoped_release>())
        .def("dump_params_json", &GCMCCPUSimulation::dumpParamsJson,
             py::arg("filename"),
             py::call_guard<py::gil_scoped_release>())
        .def("get_config", &GCMCCPUSimulation::getConfig)
        .def("update_config", &GCMCCPUSimulation::updateConfig,
             py::arg("config"))
        .def("get_statistics", [](const GCMCCPUSimulation& self) {
            return cpu_stats_to_dict(self.getStatistics());
        })
        .def("get_counters", [](const GCMCCPUSimulation& self) {
            return cpu_counters_to_dict(self.getCounters());
        })
        .def("get_fragment_info", [](const GCMCCPUSimulation& self) {
            return cpu_fragments_to_list(self.getFragmentInfo());
        })
        .def("get_last_move", [](const GCMCCPUSimulation& self) {
            return cpu_acceptance_record_to_dict(self.getLastMove());
        })
        .def("get_moves", [](const GCMCCPUSimulation& self, size_t n) {
            return cpu_acceptance_records_to_list(self.getMoves(n));
        }, py::arg("n"))
        .def("dump_lj_matrix", &GCMCCPUSimulation::dumpLJMatrix,
             py::call_guard<py::gil_scoped_release>())
        .def("print_statistics", &GCMCCPUSimulation::printStatistics)
        .def("get_acceptance_records", [](const GCMCCPUSimulation& self, size_t n) {
            return cpu_acceptance_records_to_list(self.getMoves(n));
        }, py::arg("n"));

    m.def("run_gcmc_cpu",
          &run_gcmc_cpu_config,
          py::arg("config"),
          py::arg("resumeCheckpoint") = "",
          py::arg("dumpAccept") = "",
          py::arg("dumpParams") = "",
          py::arg("diagnosticsBuffer") = 65536,
          "Run the CPU GCMC simulation through the pybind interface.");
}

} // namespace platform
} // namespace bindings
} // namespace pygcmc
