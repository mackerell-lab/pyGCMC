#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/stl_bind.h>
// Use the new MovementAPI instead of including individual headers
#include "../../platform/cpu/movement/MovementAPI.hpp"
#include "../../platform/cpu/movement/core/MovementMain.hpp"
#include "../../platform/cpu/movement/common/MovementParams.hpp"
#include "../../platform/cpu/movement/common/MovementResult.hpp"
#include "../../platform/cpu/movement/common/MovementStatistics.hpp"
#include "../../platform/cpu/movement/pool/ActivePool.hpp"
#include "../../platform/cpu/movement/common/MovementUtils.hpp"  // For Vector3
#include "../../model/montecarlo/MCMain.hpp"  // For MCState

namespace py = pybind11;
using namespace pygcmc::platform::cpu::movement;

// Forward declaration for FragmentReservoir bindings
namespace pygcmc {
namespace bindings {
namespace movement {
    void init_fragment_reservoir_bindings(pybind11::module& m);
}
}
}

namespace pygcmc {
namespace bindings {
namespace simulation {

void init_movement_bindings(py::module& m) {
    // Create movement submodule
    auto movement = m.def_submodule("movement", "GCMC Movement operations");
    
    // Bind MovementParams class - comprehensive binding of all members
    py::class_<MovementParams>(movement, "MovementParams")
        .def(py::init<>())
        .def(py::init<double>(), py::arg("temperature"))
        .def("__repr__", [](const MovementParams& p) {
            return std::string("MovementParams(T=") + std::to_string(p.temperature) + ")";
        })
        // Basic thermodynamic parameters
        .def_readwrite("temperature", &MovementParams::temperature)
        .def_readwrite("beta", &MovementParams::beta)
        .def_readwrite("chemicalPotential", &MovementParams::chemicalPotential)
        // Cavity bias parameters
        .def_readwrite("useCavityBias", &MovementParams::useCavityBias)
        .def_readwrite("cavityGridSpacing", &MovementParams::cavityGridSpacing)
        .def_readwrite("probeRadius", &MovementParams::probeRadius)
        .def_readwrite("cavityUpdateFrequency", &MovementParams::cavityUpdateFrequency)
        // Configurational bias parameters
        .def_readwrite("useConfigBias", &MovementParams::useConfigBias)
        .def_readwrite("useConfigBiasForInsertion", &MovementParams::useConfigBiasForInsertion)
        .def_readwrite("numConfigTrials", &MovementParams::numConfigTrials)
        .def_readwrite("includeTranslationInConfig", &MovementParams::includeTranslationInConfig)
        .def_readwrite("configTranslationRange", &MovementParams::configTranslationRange)
        // Movement limits
        .def_readwrite("maxTranslation", &MovementParams::maxTranslation)
        .def_readwrite("maxRotation", &MovementParams::maxRotation)
        // Numerical stability parameters
        .def_readwrite("useLogSpace", &MovementParams::useLogSpace)
        .def_readwrite("logSpaceMin", &MovementParams::logSpaceMin)
        .def_readwrite("logSpaceMax", &MovementParams::logSpaceMax)
        // Two-stage strategy parameters
        .def_readwrite("useTwoStageStrategy", &MovementParams::useTwoStageStrategy)
        .def_readwrite("pgpThreshold", &MovementParams::pgpThreshold)
        .def_readwrite("maxCandidatesStage2", &MovementParams::maxCandidatesStage2)
        // Active pool parameters
        .def_readwrite("maxAtoms", &MovementParams::maxAtoms)
        .def_readwrite("maxResidues", &MovementParams::maxResidues)
        .def_readwrite("fragmentationThreshold", &MovementParams::fragmentationThreshold)
        // System parameters
        .def_readwrite("volumeNm3", &MovementParams::volumeNm3)
        .def_readwrite("idealGasConcentration", &MovementParams::idealGasConcentration)
        .def_readwrite("thermalLambdaNm", &MovementParams::thermalLambdaNm)
        // Random number generator seed
        .def_readwrite("seed", &MovementParams::seed)
        // Movement probabilities
        .def_readwrite("insertionProbability", &MovementParams::insertionProbability)
        .def_readwrite("deletionProbability", &MovementParams::deletionProbability)
        .def_readwrite("translationProbability", &MovementParams::translationProbability)
        .def_readwrite("rotationProbability", &MovementParams::rotationProbability)
        // Multi-insertion CBMC parameters
        .def_readwrite("useMultiInsertionCBMC", &MovementParams::useMultiInsertionCBMC)
        .def_readwrite("maxParallelInsertions", &MovementParams::maxParallelInsertions)
        .def_readwrite("minRegionSeparationNm", &MovementParams::minRegionSeparationNm)
        .def_readwrite("multiDisplacementFraction", &MovementParams::multiDisplacementFraction)
        .def_readwrite("multiUseRegionVolume", &MovementParams::multiUseRegionVolume)
        // Proposal layer parameters
        .def_readwrite("proposalMode", &MovementParams::proposalMode)
        // Adaptive mode thresholds
        .def_readwrite("autoOccupancySparse", &MovementParams::autoOccupancySparse)
        .def_readwrite("autoOccupancyDense", &MovementParams::autoOccupancyDense)
        .def_readwrite("autoNcavMin", &MovementParams::autoNcavMin)
        .def_readwrite("autoFindCavMaxMs", &MovementParams::autoFindCavMaxMs)
        // Performance features
        .def_readwrite("useIncrementalCavityUpdate", &MovementParams::useIncrementalCavityUpdate)
        .def_readwrite("useStencilOptimization", &MovementParams::useStencilOptimization)
        .def_readwrite("useColorClassFastPath", &MovementParams::useColorClassFastPath)
        // Diagnostic options
        .def_readwrite("fillProposalInfo", &MovementParams::fillProposalInfo)
        // Methods
        .def("updateDerivedParameters", &MovementParams::updateDerivedParameters)
        .def("validateParameters", &MovementParams::validateParameters);
    
    // Bind MovementResult class
    py::class_<MovementResult>(movement, "MovementResult")
        .def(py::init<>())
        .def(py::init<bool, double, double, const std::string&>(),
             py::arg("accepted"), py::arg("energyChange"), 
             py::arg("acceptanceProbability"), py::arg("moveType"))
        .def_readonly("accepted", &MovementResult::accepted)
        .def_readonly("energyChange", &MovementResult::energyChange)
        .def_readonly("acceptanceProbability", &MovementResult::acceptanceProbability)
        .def_readonly("moveType", &MovementResult::moveType)
        .def_readonly("residueIndex", &MovementResult::residueIndex)
        .def_readonly("moleculeType", &MovementResult::moleculeType)
        .def_readonly("cavityBiasFactor", &MovementResult::cavityBiasFactor)
        .def_readonly("configBiasFactor", &MovementResult::configBiasFactor)
        .def_readonly("numConfigTrials", &MovementResult::numConfigTrials)
        .def_readonly("computeTimeMs", &MovementResult::computeTimeMs)
        .def_readonly("rejectReason", &MovementResult::rejectReason)
        .def_readonly("numericalError", &MovementResult::numericalError)
        .def_readonly("proposalInfoFilled", &MovementResult::proposalInfoFilled)
        // 实现为方法以保持向后兼容性（原simulation模块中是方法）
        .def("isSuccessful", [](const MovementResult& r) { return r.accepted; }, 
             "Check if the move was successful")
        .def("summary", [](const MovementResult& r) {
            return std::string(r.moveType) + " " + (r.accepted ? "accepted" : "rejected") + 
                   " (ΔE=" + std::to_string(r.energyChange) + " kJ/mol)";
        }, "Get a summary string of the result")
        .def("__repr__", [](const MovementResult& r) {
            return std::string("MovementResult(") + 
                   (r.accepted ? "accepted" : "rejected") + 
                   ", energy=" + std::to_string(r.energyChange) + ")";
        });
    
    // Bind MovementModule class
    py::class_<MovementModule>(movement, "MovementModule")
        .def(py::init<>())
        .def(py::init<const MovementParams&>())
        .def("setParams", &MovementModule::setParams)
        .def("getParams", &MovementModule::getParams)
        .def("attemptInsertion", &MovementModule::attemptInsertion,
             py::arg("state"), py::arg("moleculeType") = 0)
        .def("attemptDeletion", &MovementModule::attemptDeletion,
             py::arg("state"), py::arg("residueIndex") = -1)
        .def("attemptTranslation", &MovementModule::attemptTranslation,
             py::arg("state"), py::arg("residueIndex") = -1)
        .def("attemptRotation", &MovementModule::attemptRotation,
             py::arg("state"), py::arg("residueIndex") = -1)
        .def("attemptCavityBiasInsertion", &MovementModule::attemptCavityBiasInsertion,
             py::arg("state"), py::arg("moleculeType") = 0)
        .def("attemptConfigBiasRotation", &MovementModule::attemptConfigBiasRotation,
             py::arg("state"), py::arg("residueIndex") = -1)
        .def("findCavities", &MovementModule::findCavities)
        .def("calculateCavityVolume", &MovementModule::calculateCavityVolume)
        .def("calculateAcceptanceRate", &MovementModule::calculateAcceptanceRate)
        .def("getStatistics", &MovementModule::getStatistics,
             py::return_value_policy::copy)
        .def("resetStatistics", &MovementModule::resetStatistics)
        .def("attemptMultiInsertionCBMC", &MovementModule::attemptMultiInsertionCBMC,
             py::arg("state"), py::arg("moleculeType") = 0)
        .def("getProposalStats", &MovementModule::getProposalStatsMap,
             "Get proposal statistics as a map")
        .def("getCavityStats", &MovementModule::getCavityStatsMap,
             "Get cavity statistics as a map")
        .def("__repr__", [](const MovementModule&) {
            return std::string("MovementModule(temperature=300.0)");
        });
    
    // Basic GCMC functions using MovementAPI
    movement.def("initializeGCMC", &initializeGCMC,
                "Initialize GCMC module with state");
    
    movement.def("performGCMCMove", &performGCMCMove,
                "Perform a single GCMC move");
    
    movement.def("runGCMCSteps", &runGCMCSteps,
                "Run multiple GCMC steps",
                py::arg("nSteps"));
    
    movement.def("getGCMCStatistics", &getGCMCStatistics,
                "Get GCMC statistics");
    
    // Fragment reservoir functions
    movement.def("addFragmentToReservoir", &addFragmentToReservoir,
                "Add a fragment to the reservoir",
                py::arg("name"),
                py::arg("atoms"),
                py::arg("chemicalPotential") = -15.7);
    
    movement.def("clearFragmentReservoir", &clearFragmentReservoir,
                "Clear the fragment reservoir");
    
    movement.def("getFragmentCount", &getFragmentCount,
                "Get number of fragments in reservoir");
    
    // Individual move types
    movement.def("attemptInsertion", &attemptInsertion,
                "Attempt an insertion move");
    
    movement.def("attemptDeletion", &attemptDeletion,
                "Attempt a deletion move");
    
    movement.def("attemptTranslation", &attemptTranslation,
                "Attempt a translation move",
                py::arg("state"),
                py::arg("residueIndex"));
    
    movement.def("attemptRotation", &attemptRotation,
                "Attempt a rotation move",
                py::arg("state"),
                py::arg("residueIndex"));
    
    // Utility functions
    movement.def("setRandomSeed", &setRandomSeed,
                "Set random seed for movement operations",
                py::arg("seed"));
    
    movement.def("setCavityBiasEnabled", &setCavityBiasEnabled,
                "Enable/disable cavity bias",
                py::arg("enabled"));
    
    movement.def("setConfigBiasEnabled", &setConfigBiasEnabled,
                "Enable/disable configurational bias",
                py::arg("enabled"));
    
    // Bind Vector3 class (required for several other classes)
    // Note: using the movement namespace Vector3
    py::class_<pygcmc::platform::cpu::movement::Vector3>(movement, "Vector3")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def_readwrite("x", &pygcmc::platform::cpu::movement::Vector3::x)
        .def_readwrite("y", &pygcmc::platform::cpu::movement::Vector3::y)
        .def_readwrite("z", &pygcmc::platform::cpu::movement::Vector3::z);
    
    // Bind ActivePool::ResidueMetadata
    py::class_<ActivePool::ResidueMetadata>(movement, "ResidueMetadata")
        .def(py::init<>())
        .def_readwrite("residueType", &ActivePool::ResidueMetadata::residueType)
        .def_readwrite("active", &ActivePool::ResidueMetadata::active)
        .def_readwrite("atomStartIndex", &ActivePool::ResidueMetadata::atomStartIndex)
        .def_readwrite("atomCount", &ActivePool::ResidueMetadata::atomCount)
        .def_readwrite("centerOfMass", &ActivePool::ResidueMetadata::centerOfMass)
        .def_readwrite("insertionTime", &ActivePool::ResidueMetadata::insertionTime);
    
    // Bind ActivePool::Statistics
    py::class_<ActivePool::Statistics>(movement, "ActivePoolStatistics")
        .def(py::init<>())
        .def_readonly("totalInserts", &ActivePool::Statistics::totalInserts)
        .def_readonly("totalDeletes", &ActivePool::Statistics::totalDeletes)
        .def_readonly("compactions", &ActivePool::Statistics::compactions)
        .def_readonly("peakAtoms", &ActivePool::Statistics::peakAtoms)
        .def_readonly("peakResidues", &ActivePool::Statistics::peakResidues)
        .def_readonly("averageFragmentation", &ActivePool::Statistics::averageFragmentation)
        .def_readonly("batchOperations", &ActivePool::Statistics::batchOperations);
    
    // Bind ActivePool class
    py::class_<ActivePool>(movement, "ActivePool")
        .def(py::init<int, int>(), py::arg("maxAtoms") = 100000, py::arg("maxResidues") = 30000)
        .def("insertMolecule", &ActivePool::insertMolecule,
             py::arg("atoms"), py::arg("resType") = 0)
        .def("deleteResidue", &ActivePool::deleteResidue)
        .def("compact", &ActivePool::compact, py::arg("force") = false)
        .def("syncToState", &ActivePool::syncToState)
        .def("syncFromState", &ActivePool::syncFromState)
        .def("getFragmentation", &ActivePool::getFragmentation)
        .def("getActiveCounts", &ActivePool::getActiveCounts)
        .def("getActiveResidueIndices", &ActivePool::getActiveResidueIndices)
        .def("isResidueActive", &ActivePool::isResidueActive)
        .def("getStatistics", &ActivePool::getStatistics,
             py::return_value_policy::reference_internal)
        .def("resetStatistics", &ActivePool::resetStatistics)
        .def("setFragmentationThreshold", &ActivePool::setFragmentationThreshold)
        .def("getFragmentationThreshold", &ActivePool::getFragmentationThreshold)
        .def("getMaxAtoms", &ActivePool::getMaxAtoms)
        .def("getMaxResidues", &ActivePool::getMaxResidues)
        .def("canInsert", &ActivePool::canInsert);
    
    // Bind MovementStatistics class 
    py::class_<MovementStatistics>(movement, "MovementStatistics")
        .def(py::init<>())
        .def_readwrite("attempts", &MovementStatistics::attempts)
        .def_readwrite("accepts", &MovementStatistics::accepts)
        .def_readwrite("totalEnergyChange", &MovementStatistics::totalEnergyChange)
        .def("acceptanceRate", &MovementStatistics::acceptanceRate);
    
    // Register FragmentReservoir and related classes
    pygcmc::bindings::movement::init_fragment_reservoir_bindings(movement);
}

} // namespace simulation
} // namespace bindings
} // namespace pygcmc