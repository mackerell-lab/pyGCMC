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
#include "../../platform/cpu/movement/bias/CavityBiasCore.hpp"  // For CavityBiasCore

namespace py = pybind11;
using namespace ::pygcmc::platform::cpu::movement;

// Forward declaration for FragmentReservoir bindings
namespace pygcmc {
namespace bindings {
namespace platform {
    void init_fragment_reservoir_bindings(pybind11::module& m);
}
}
}

namespace pygcmc {
namespace bindings {
namespace platform {

void init_movement_bindings(py::module& m) {
    // Create movement submodule
    auto movement = m.def_submodule("movement", "GCMC Movement operations");
    
    // Expose unified grand-canonical bookkeeping types
    py::class_<gcmc::GrandCanonicalTerms>(movement, "GrandCanonicalTerms")
        .def(py::init<>())
        .def_readwrite("speciesId", &gcmc::GrandCanonicalTerms::speciesId)
        .def_readwrite("countBefore", &gcmc::GrandCanonicalTerms::countBefore)
        .def_readwrite("countAfter", &gcmc::GrandCanonicalTerms::countAfter)
        .def_readwrite("beta", &gcmc::GrandCanonicalTerms::beta)
        .def_readwrite("chemicalPotential", &gcmc::GrandCanonicalTerms::chemicalPotential)
        .def_readwrite("deltaEnergy", &gcmc::GrandCanonicalTerms::deltaEnergy)
        .def_readwrite("logVolume", &gcmc::GrandCanonicalTerms::logVolume)
        .def_readwrite("logLambda3", &gcmc::GrandCanonicalTerms::logLambda3)
        .def_readwrite("logProposalForward", &gcmc::GrandCanonicalTerms::logProposalForward)
        .def_readwrite("logProposalReverse", &gcmc::GrandCanonicalTerms::logProposalReverse)
        .def_readwrite("logCavityForward", &gcmc::GrandCanonicalTerms::logCavityForward)
        .def_readwrite("logCavityReverse", &gcmc::GrandCanonicalTerms::logCavityReverse)
        .def_readwrite("logRosenbluthForward", &gcmc::GrandCanonicalTerms::logRosenbluthForward)
        .def_readwrite("logRosenbluthReverse", &gcmc::GrandCanonicalTerms::logRosenbluthReverse)
        .def_readwrite("logExtraForward", &gcmc::GrandCanonicalTerms::logExtraForward)
        .def_readwrite("logExtraReverse", &gcmc::GrandCanonicalTerms::logExtraReverse);

    py::class_<gcmc::GrandCanonicalEvaluation>(movement, "GrandCanonicalEvaluation")
        .def(py::init<>())
        .def_readwrite("probability", &gcmc::GrandCanonicalEvaluation::probability)
        .def_readwrite("logRatio", &gcmc::GrandCanonicalEvaluation::logRatio);

    py::class_<MovementParams::CavityFragmentConfig>(movement, "CavityFragmentConfig")
        .def(py::init<>())
        .def_readwrite("gridSpacingNm", &MovementParams::CavityFragmentConfig::gridSpacingNm)
        .def_readwrite("probeRadiusNm", &MovementParams::CavityFragmentConfig::probeRadiusNm)
        .def_readwrite("maskId", &MovementParams::CavityFragmentConfig::maskId);

    py::class_<MovementParams::MoveProbabilitySet>(movement, "MoveProbabilitySet")
        .def(py::init<>())
        .def_readwrite("insertion", &MovementParams::MoveProbabilitySet::insertion)
        .def_readwrite("deletion", &MovementParams::MoveProbabilitySet::deletion)
        .def_readwrite("translation", &MovementParams::MoveProbabilitySet::translation)
        .def_readwrite("rotation", &MovementParams::MoveProbabilitySet::rotation);

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
        .def_readwrite("attemptProbInsertion", &MovementParams::attemptProbInsertion)
        .def_readwrite("attemptProbDeletion", &MovementParams::attemptProbDeletion)
        .def_readwrite("attemptProbTranslation", &MovementParams::attemptProbTranslation)
        .def_readwrite("attemptProbRotation", &MovementParams::attemptProbRotation)
        .def_readwrite("targetFragmentCounts", &MovementParams::targetFragmentCounts)
        .def_readwrite("removeInitFlags", &MovementParams::removeInitFlags)
        .def_readwrite("removeExcessFlags", &MovementParams::removeExcessFlags)
        .def_readwrite("excessRemovalThresholds", &MovementParams::excessRemovalThresholds)
        .def_readwrite("defaultExcessThreshold", &MovementParams::defaultExcessThreshold)
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
        .def_readwrite("cavityFragmentConfigs", &MovementParams::cavityFragmentConfigs)
        // Methods
        .def("get_move_probability_set",
             &MovementParams::getMoveProbabilitySet,
             py::arg("fragment_type"))
        .def("get_move_probability_set_biased",
             &MovementParams::getBiasedMoveProbabilitySet,
             py::arg("fragment_type"),
             py::arg("population_before"))
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
        .def_readonly("cbmcTrialsUsed", &MovementResult::cbmcTrialsUsed)
        .def_readonly("logAcceptanceRatio", &MovementResult::logAcceptanceRatio)
        .def_readonly("logProposalForward", &MovementResult::logProposalForward)
        .def_readonly("logProposalReverse", &MovementResult::logProposalReverse)
        .def_readonly("volumeNm3", &MovementResult::volumeNm3)
        .def_readonly("effectiveVolumeNm3", &MovementResult::effectiveVolumeNm3)
        .def_readonly("cavityVolumeNm3", &MovementResult::cavityVolumeNm3)
        .def_readonly("logVolume", &MovementResult::logVolume)
        .def_readonly("logCavityFactor", &MovementResult::logCavityFactor)
        .def_readonly("lambdaNm", &MovementResult::lambdaNm)
        .def_readonly("logLambda3", &MovementResult::logLambda3)
        .def_readonly("logWForward", &MovementResult::logWForward)
        .def_readonly("logWReverse", &MovementResult::logWReverse)
        .def_readonly("hasGrandTerms", &MovementResult::hasGrandTerms)
        .def_readonly("grandTerms", &MovementResult::grandTerms)
        .def_readonly("grandEvaluation", &MovementResult::grandEvaluation)
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
        .def("getPopulationControlStatsMap", &MovementModule::getPopulationControlStatsMap,
             "Get population control enforcement counters")
        .def("__repr__", [](const MovementModule&) {
            return std::string("MovementModule(temperature=300.0)");
        });

    // Expose acceptance probability helpers for regression tests
    movement.def(
        "calculate_insertion_probability_basic",
        [](int n, double deltaE, double beta, double chemPotential,
           double cavityBias, double volumeNm3) {
            return utils::LogSpaceCalculator::calculateInsertionProbability(
                n, deltaE, beta, chemPotential, cavityBias, volumeNm3, true);
        },
        py::arg("n"),
        py::arg("deltaE"),
        py::arg("beta"),
        py::arg("chemicalPotential"),
        py::arg("cavityBias"),
        py::arg("volumeNm3"),
        "Compute insertion acceptance probability (log-space) including cavity bias");

    movement.def(
        "calculate_insertion_probability_with_lambda",
        [](int n, double deltaE, double beta, double chemPotential,
           double cavityBias, double volumeNm3, double thermalLambdaNm) {
            return utils::LogSpaceCalculator::calculateInsertionProbabilityWithLambda(
                n, deltaE, beta, chemPotential, cavityBias, volumeNm3,
                thermalLambdaNm, true);
        },
        py::arg("n"),
        py::arg("deltaE"),
        py::arg("beta"),
        py::arg("chemicalPotential"),
        py::arg("cavityBias"),
        py::arg("volumeNm3"),
        py::arg("thermalLambdaNm"),
        "Compute insertion acceptance probability with thermal wavelength factor");

    movement.def(
        "calculate_deletion_probability_basic",
        [](int n, double deltaE, double beta, double chemPotential,
           double volumeNm3) {
            return utils::LogSpaceCalculator::calculateDeletionProbability(
                n, deltaE, beta, chemPotential, volumeNm3, true);
        },
        py::arg("n"),
        py::arg("deltaE"),
        py::arg("beta"),
        py::arg("chemicalPotential"),
        py::arg("volumeNm3"),
        "Compute deletion acceptance probability without cavity or lambda terms");

    movement.def(
        "calculate_deletion_probability_with_cavity",
        [](int n, double deltaE, double beta, double chemPotential,
           double cavityBias, double volumeNm3) {
            return utils::LogSpaceCalculator::calculateDeletionProbabilityWithCavity(
                n, deltaE, beta, chemPotential, cavityBias, volumeNm3, true);
        },
        py::arg("n"),
        py::arg("deltaE"),
        py::arg("beta"),
        py::arg("chemicalPotential"),
        py::arg("cavityBias"),
        py::arg("volumeNm3"),
        "Compute deletion acceptance probability including cavity bias");

    movement.def(
        "calculate_deletion_probability_with_cavity_and_lambda",
        [](int n, double deltaE, double beta, double chemPotential,
           double cavityBias, double volumeNm3, double thermalLambdaNm) {
            return utils::LogSpaceCalculator::calculateDeletionProbabilityWithCavityAndLambda(
                n, deltaE, beta, chemPotential, cavityBias, volumeNm3,
                thermalLambdaNm, true);
        },
        py::arg("n"),
        py::arg("deltaE"),
        py::arg("beta"),
        py::arg("chemicalPotential"),
        py::arg("cavityBias"),
        py::arg("volumeNm3"),
        py::arg("thermalLambdaNm"),
        "Compute deletion acceptance probability with both cavity and thermal wavelength factors");
    
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
    py::class_<::pygcmc::platform::cpu::movement::Vector3>(movement, "Vector3")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def_readwrite("x", &::pygcmc::platform::cpu::movement::Vector3::x)
        .def_readwrite("y", &::pygcmc::platform::cpu::movement::Vector3::y)
        .def_readwrite("z", &::pygcmc::platform::cpu::movement::Vector3::z);
    
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
    
    // Bind CavityMode enum
    py::enum_<CavityMode>(movement, "CavityMode")
        .value("FAST_APPROX", CavityMode::FAST_APPROX, "Global cavity fraction (default)")
        .value("CLUSTER_VOLUME", CavityMode::CLUSTER_VOLUME, "Cluster-based accurate sampling")
        .value("LOCAL_VEFF", CavityMode::LOCAL_VEFF, "Local effective volume (highest accuracy)")
        .export_values();
    
    // Bind CavityBiasCore class
    py::class_<CavityBiasCore>(movement, "CavityBiasCore")
        .def(py::init<double, double>(), 
             py::arg("gridSpacing") = 0.25, 
             py::arg("probeRadius") = 0.14,
             "Create CavityBiasCore with grid spacing and probe radius in nm")
        .def("calculateCavityVolume",
             [](CavityBiasCore& core,
                const model::montecarlo::MCState& state,
                CavityMode mode,
                int speciesId) {
                 return core.calculateCavityVolume(state, mode, speciesId);
             },
             py::arg("state"),
             py::arg("mode"),
             py::arg("speciesId") = -1,
             "Calculate cavity volume in nm^3")
        .def("proposeCavityPosition", &CavityBiasCore::proposeCavityPosition,
             py::arg("state"),
             py::arg("mode"),
             "Propose a position inside a cavity")
        .def("invalidateCache", &CavityBiasCore::invalidateCache,
             "Invalidate the internal cache")
        .def("setGridSpacing", &CavityBiasCore::setGridSpacing,
             py::arg("spacing"),
             "Set grid spacing in nm")
        .def("setProbeRadius", &CavityBiasCore::setProbeRadius,
             py::arg("radius"),
             "Set probe radius in nm")
        .def("setSpeciesParameters", &CavityBiasCore::setSpeciesParameters,
             py::arg("speciesId"),
             py::arg("gridSpacing"),
             py::arg("probeRadius"),
             py::arg("maskId") = -1,
             "Override parameters for a specific species index")
        .def("clearSpeciesParameters", &CavityBiasCore::clearSpeciesParameters,
             "Clear species-specific overrides")
        .def("__repr__", [](const CavityBiasCore& /*c*/) {
            return std::string("CavityBiasCore(gridSpacing=0.25, probeRadius=0.14)");
        });
    
    // Register FragmentReservoir and related classes
    pygcmc::bindings::platform::init_fragment_reservoir_bindings(movement);
}

} // namespace platform
} // namespace bindings
} // namespace pygcmc
