#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "platform/cpu/movement/MovementModule.hpp"
#include "model/montecarlo/MCMain.hpp"  // For MCState
#include "platform/cpu/movement/reservoir/fragment_reservoir.hpp"  // For FragmentReservoir
// All necessary types are included via MovementModule.hpp

namespace py = pybind11;
using namespace pygcmc::platform::cpu::movement;

namespace pygcmc {
namespace bindings {

// Forward declaration for fragment reservoir bindings
namespace movement {
    void init_fragment_reservoir_bindings(py::module& m);
}

namespace simulation {

void init_movement_bindings(py::module& m) {
    // Create movement submodule
    auto movement = m.def_submodule("movement", "GCMC Movement operations");
    
    // Bind Vector3 class
    py::class_<Vector3>(movement, "Vector3")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def_readwrite("x", &Vector3::x)
        .def_readwrite("y", &Vector3::y)
        .def_readwrite("z", &Vector3::z)
        .def("__repr__", [](const Vector3& v) {
            return "<Vector3(" + std::to_string(v.x) + ", " + 
                   std::to_string(v.y) + ", " + std::to_string(v.z) + ")>";
        });
    
    // Bind MovementParams
    py::class_<MovementParams>(movement, "MovementParams")
        .def(py::init<>())
        .def(py::init<double>(), py::arg("temperature"))
        .def_readwrite("temperature", &MovementParams::temperature, 
                      "Temperature in Kelvin")
        .def_readwrite("beta", &MovementParams::beta, 
                      "1/kT in mol/kJ")
        .def_readwrite("chemicalPotential", &MovementParams::chemicalPotential,
                      "Chemical potential in kJ/mol")
        .def_readwrite("useCavityBias", &MovementParams::useCavityBias,
                      "Enable cavity bias for insertion")
        .def_readwrite("cavityGridSpacing", &MovementParams::cavityGridSpacing,
                      "Grid spacing for cavity detection in nm")
        .def_readwrite("probeRadius", &MovementParams::probeRadius,
                      "Probe radius for cavity detection in nm")
        .def_readwrite("useConfigBias", &MovementParams::useConfigBias,
                      "Enable configurational bias")
        .def_readwrite("useConfigBiasForInsertion", &MovementParams::useConfigBiasForInsertion,
                      "Enable CBMC for insertion moves (two-step method)")
        .def_readwrite("numConfigTrials", &MovementParams::numConfigTrials,
                      "Number of trial configurations")
        .def_readwrite("configTranslationRange", &MovementParams::configTranslationRange,
                      "Translation range for trial configurations in nm")
        .def_readwrite("maxTranslation", &MovementParams::maxTranslation,
                      "Maximum translation distance in nm")
        .def_readwrite("maxRotation", &MovementParams::maxRotation,
                      "Maximum rotation angle in radians")
        .def_readwrite("useLogSpace", &MovementParams::useLogSpace,
                      "Use log-space calculations for numerical stability")
        .def_readwrite("maxAtoms", &MovementParams::maxAtoms,
                      "Maximum atoms in active pool")
        .def_readwrite("maxResidues", &MovementParams::maxResidues,
                      "Maximum residues in active pool")
        .def_readwrite("seed", &MovementParams::seed,
                      "Random number generator seed (0 = use time-based seed)")
        .def_readwrite("thermalLambdaNm", &MovementParams::thermalLambdaNm,
                      "Thermal de Broglie wavelength in nm (default 1.0 for compatibility)")
        .def_readwrite("useMultiInsertionCBMC", &MovementParams::useMultiInsertionCBMC,
                      "Enable multi-insertion CBMC")
        .def_readwrite("maxParallelInsertions", &MovementParams::maxParallelInsertions,
                      "Maximum number of parallel insertions")
        .def_readwrite("minRegionSeparationNm", &MovementParams::minRegionSeparationNm,
                      "Minimum separation between insertion regions in nm")
        .def_readwrite("multiDisplacementFraction", &MovementParams::multiDisplacementFraction,
                      "Sampling radius fraction in region (0..1]")
        .def_readwrite("multiUseRegionVolume", &MovementParams::multiUseRegionVolume,
                      "Use region volume for Veff when true; otherwise box volume")
        .def_readwrite("proposalMode", &MovementParams::proposalMode,
                      "Proposal sampling mode: 0=Uniform, 1=Cavity, 2=Color, 3=Cluster, 4=Adaptive "
                      "(requires PYGCMC_USE_PROPOSAL_LAYER compile flag)")
        // P2: Adaptive mode thresholds
        .def_readwrite("autoOccupancySparse", &MovementParams::autoOccupancySparse,
                      "Switch to Uniform below this occupancy (default 0.3)")
        .def_readwrite("autoOccupancyDense", &MovementParams::autoOccupancyDense,
                      "Switch to Cluster above this occupancy (default 0.7)")
        .def_readwrite("autoNcavMin", &MovementParams::autoNcavMin,
                      "Minimum cavities for Cavity mode (default 100)")
        .def_readwrite("autoFindCavMaxMs", &MovementParams::autoFindCavMaxMs,
                      "Max time (ms) before switching modes (default 10.0)")
        // P3: Optional performance features
        .def_readwrite("useIncrementalCavityUpdate", &MovementParams::useIncrementalCavityUpdate,
                      "Enable incremental cavity update on accept (default false)")
        .def_readwrite("useStencilOptimization", &MovementParams::useStencilOptimization,
                      "Use precomputed sphere stencils (default true)")
        .def_readwrite("useColorClassFastPath", &MovementParams::useColorClassFastPath,
                      "Use color class index tables (default false)")
        .def_readwrite("fillProposalInfo", &MovementParams::fillProposalInfo,
                      "Fill detailed proposal info in result (default false)")
        .def("updateDerivedParameters", &MovementParams::updateDerivedParameters,
                      "Update derived parameters after changing temperature")
        .def("__repr__", [](const MovementParams& p) {
            return "<MovementParams T=" + std::to_string(p.temperature) + 
                   " μ=" + std::to_string(p.chemicalPotential) + 
                   " multiInsert=" + (p.useMultiInsertionCBMC ? "ON" : "OFF") + ">";
        });
    
    // Bind MovementResult
    py::class_<MovementResult>(movement, "MovementResult")
        .def(py::init<>())
        .def(py::init<bool, double, double, const std::string&>(),
             py::arg("accepted"), py::arg("energyChange"), 
             py::arg("acceptanceProbability"), py::arg("moveType"))
        .def_readonly("accepted", &MovementResult::accepted,
                     "Whether the move was accepted")
        .def_readonly("energyChange", &MovementResult::energyChange,
                     "Change in energy in kJ/mol")
        .def_readonly("acceptanceProbability", &MovementResult::acceptanceProbability,
                     "Calculated acceptance probability")
        .def_readonly("moveType", &MovementResult::moveType,
                     "Type of move performed")
        .def_readonly("residueIndex", &MovementResult::residueIndex,
                     "Index of affected residue")
        .def_readonly("cavityBiasFactor", &MovementResult::cavityBiasFactor,
                     "Cavity bias correction factor")
        .def_readonly("configBiasFactor", &MovementResult::configBiasFactor,
                     "Configurational bias correction factor")
        .def_readonly("usedCavity", &MovementResult::usedCavity,
                     "Whether cavity was used in proposal")
        .def_readonly("mproposal", &MovementResult::mproposal,
                     "Number of proposal positions (-1 if not applicable)")
        .def_readonly("vregion", &MovementResult::vregion,
                     "Region volume in nm³ (-1 if not applicable)")
        // P2: Enhanced diagnostics
        .def_readonly("proposalMode", &MovementResult::proposalMode,
                     "Current proposal mode (0-4, -1 if not set)")
        .def_readonly("selectedType", &MovementResult::selectedType,
                     "ProposalType enum value actually used")
        .def_readonly("proposalTimeMs", &MovementResult::proposalTimeMs,
                     "Time for proposal generation (ms)")
        .def_readonly("findCavTimeMs", &MovementResult::findCavTimeMs,
                     "Time for cavity finding (ms)")
        // Latest proposal info
        .def_readonly("proposalNorm", &MovementResult::proposalNorm,
                     "q_norm for detailed balance diagnostics")
        .def_readonly("proposalPosX", &MovementResult::proposalPosX,
                     "Proposed position X (nm)")
        .def_readonly("proposalPosY", &MovementResult::proposalPosY,
                     "Proposed position Y (nm)")
        .def_readonly("proposalPosZ", &MovementResult::proposalPosZ,
                     "Proposed position Z (nm)")
        .def_readonly("proposalInfoFilled", &MovementResult::proposalInfoFilled,
                     "Whether proposal info was filled")
        .def_readonly("computeTimeMs", &MovementResult::computeTimeMs,
                     "Time taken for the move in milliseconds")
        .def("isSuccessful", &MovementResult::isSuccessful,
                     "Check if move was successful")
        .def("summary", &MovementResult::summary,
                     "Get a summary string of the result")
        .def("__repr__", [](const MovementResult& r) {
            return "<MovementResult " + r.moveType + ": " + 
                   (r.accepted ? "accepted" : "rejected") + ">";
        });
    
    // Bind MovementModule
    py::class_<MovementModule>(movement, "MovementModule", 
                               "Main GCMC Movement Module for CPU platform")
        .def(py::init<>())
        .def(py::init<const MovementParams&>(), py::arg("params"))
        
        // Main movement functions
        .def("attemptInsertion", &MovementModule::attemptInsertion,
             py::arg("state"), py::arg("moleculeType") = 0,
             "Attempt to insert a molecule into the system")
        .def("attemptDeletion", &MovementModule::attemptDeletion,
             py::arg("state"), py::arg("residueIndex") = -1,
             "Attempt to delete a molecule from the system")
        .def("attemptTranslation", &MovementModule::attemptTranslation,
             py::arg("state"), py::arg("residueIndex") = -1,
             "Attempt to translate a molecule")
        .def("attemptRotation", &MovementModule::attemptRotation,
             py::arg("state"), py::arg("residueIndex") = -1,
             "Attempt to rotate a molecule")
        
        // Biased moves
        .def("attemptCavityBiasInsertion", &MovementModule::attemptCavityBiasInsertion,
             py::arg("state"), py::arg("moleculeType") = 0,
             "Attempt insertion with cavity bias")
        .def("attemptConfigBiasRotation", &MovementModule::attemptConfigBiasRotation,
             py::arg("state"), py::arg("residueIndex") = -1,
             "Attempt rotation with configurational bias")
        .def("attemptMultiInsertionCBMC", &MovementModule::attemptMultiInsertionCBMC,
             py::arg("state"), py::arg("moleculeType") = 0,
             "Attempt multiple parallel insertions using CBMC")
        
        // Utilities
        .def("findCavities", &MovementModule::findCavities,
             py::arg("state"),
             "Find cavities in the system",
             py::return_value_policy::move)
        .def("calculateAcceptanceRate", &MovementModule::calculateAcceptanceRate,
             py::arg("moveType"),
             "Calculate acceptance rate for a move type")
        .def("resetStatistics", &MovementModule::resetStatistics,
             "Reset all statistics")
        
        // Configuration
        .def("setParams", &MovementModule::setParams,
             py::arg("params"),
             "Set movement parameters")
        .def("getParams", &MovementModule::getParams,
             "Get current movement parameters")
        
        // Statistics access
        .def("getStatistics", &MovementModule::getStatistics,
             py::return_value_policy::reference_internal,
             "Get movement statistics")
        
        // P2: Enhanced statistics access with structured export
        .def("getProposalStats", [](const MovementModule& m) {
#ifdef PYGCMC_USE_PROPOSAL_LAYER
            // When proposal layer is enabled, return structured dict
            auto statsMap = m.getProposalStatsMap();
            py::dict result;
            
            // Basic stats
            result["total_attempts"] = statsMap["total_attempts"];
            result["total_accepts"] = statsMap["total_accepts"];
            result["acceptance_rate"] = statsMap["acceptance_rate"];
            result["current_mode"] = statsMap["current_mode"];
            result["mode_transitions"] = statsMap["mode_transitions"];
            result["auto_switches"] = statsMap["auto_switches"];
            
            // Per-mode stats in nested dict
            py::dict modes;
            for (const auto& mode : {"uniform", "cavity", "color", "cluster", "adaptive"}) {
                py::dict modeStats;
                std::string prefix = std::string(mode) + "_";
                modeStats["attempts"] = statsMap[prefix + "attempts"];
                modeStats["accepts"] = statsMap[prefix + "accepts"]; 
                modeStats["accept_rate"] = statsMap[prefix + "accept_rate"];
                modes[mode] = modeStats;
            }
            result["modes"] = modes;
            
            // Timing stats in nested dict
            py::dict timings;
            py::dict proposalTime;
            proposalTime["p50"] = statsMap["proposal_time_p50_ms"];
            proposalTime["p90"] = statsMap["proposal_time_p90_ms"];
            proposalTime["p95"] = statsMap.count("proposal_time_p95_ms") ? 
                statsMap.at("proposal_time_p95_ms") : -1.0;
            proposalTime["p99"] = statsMap.count("proposal_time_p99_ms") ? 
                statsMap.at("proposal_time_p99_ms") : -1.0;
            timings["proposal_ms"] = proposalTime;
            
            py::dict findcavTime;
            findcavTime["p50"] = statsMap["findcav_time_p50_ms"];
            findcavTime["p90"] = statsMap["findcav_time_p90_ms"];
            findcavTime["p95"] = statsMap.count("findcav_time_p95_ms") ? 
                statsMap.at("findcav_time_p95_ms") : -1.0;
            findcavTime["p99"] = statsMap.count("findcav_time_p99_ms") ? 
                statsMap.at("findcav_time_p99_ms") : -1.0;
            timings["findcav_ms"] = findcavTime;
            result["timings"] = timings;
            
            // Fallback stats
            py::dict fallbacks;
            fallbacks["no_cavities"] = statsMap["fallback_no_cavities"];
            fallbacks["timeout"] = statsMap["fallback_timeout"];
            fallbacks["invalid_mode"] = statsMap["fallback_invalid_mode"];
            result["fallbacks"] = fallbacks;
            
            // Cavity info
            if (statsMap.count("cavity_count")) {
                result["cavity_count"] = statsMap["cavity_count"];
            }
            if (statsMap.count("last_ncav")) {
                result["last_ncav"] = statsMap["last_ncav"];
                result["last_occupancy"] = statsMap["last_occupancy"];
            }
            
            return result;
#else
            // Fallback to simple map conversion when proposal layer not enabled
            auto statsMap = m.getProposalStatsMap();
            py::dict result;
            for (const auto& [key, value] : statsMap) {
                result[key.c_str()] = value;
            }
            return result;
#endif
        }, "Get proposal layer statistics as structured dict (with PYGCMC_USE_PROPOSAL_LAYER) or flat dict")
        .def("getCavityStats", [](const MovementModule& m) {
            // Convert map to Python dict
            auto statsMap = m.getCavityStatsMap();
            py::dict result;
            for (const auto& [key, value] : statsMap) {
                result[key.c_str()] = value;
            }
            return result;
        }, "Get cavity manager statistics as dict")
        
        // String representation
        .def("__repr__", [](const MovementModule& m) {
            auto params = m.getParams();
            return "<MovementModule T=" + std::to_string(params.temperature) + 
                   " cavityBias=" + (params.useCavityBias ? "ON" : "OFF") +
                   " configBias=" + (params.useConfigBias ? "ON" : "OFF") + ">";
        });
    
    // Bind MovementModule::Statistics
    py::class_<MovementModule::Statistics>(movement, "MovementStatistics")
        .def_readonly("attempts", &MovementModule::Statistics::attempts)
        .def_readonly("accepts", &MovementModule::Statistics::accepts)
        .def_readonly("totalEnergyChange", &MovementModule::Statistics::totalEnergyChange)
        .def("acceptanceRate", &MovementModule::Statistics::acceptanceRate,
             "Calculate acceptance rate")
        .def_property_readonly("rate", [](const MovementModule::Statistics& s) {
            return s.acceptanceRate();
        }, "Acceptance rate as a property (alias for acceptanceRate())")
        .def("__repr__", [](const MovementModule::Statistics& s) {
            return "<Statistics attempts=" + std::to_string(s.attempts) + 
                   " accepts=" + std::to_string(s.accepts) + 
                   " rate=" + std::to_string(s.acceptanceRate()) + ">";
        });
    
    // Bind ActivePool
    py::class_<ActivePool>(movement, "ActivePool",
                          "Active pool for managing pre-allocated memory")
        .def(py::init<int, int>(), 
             py::arg("maxAtoms") = 100000, 
             py::arg("maxResidues") = 30000)
        
        // Core operations
        .def("insertMolecule", &ActivePool::insertMolecule,
             py::arg("atoms"), py::arg("resType") = 0,
             "Insert a molecule into the pool")
        .def("deleteResidue", &ActivePool::deleteResidue,
             py::arg("resIdx"),
             "Delete a residue from the pool")
        .def("compact", &ActivePool::compact,
             py::arg("force") = false,
             "Compact the pool to reduce fragmentation")
        
        // Query operations
        .def("getFragmentation", &ActivePool::getFragmentation,
             "Get current fragmentation level")
        .def("getActiveCounts", &ActivePool::getActiveCounts,
             "Get active atom and residue counts")
        .def("getActiveResidueIndices", &ActivePool::getActiveResidueIndices,
             "Get list of active residue indices")
        .def("isResidueActive", &ActivePool::isResidueActive,
             py::arg("resIdx"),
             "Check if a residue is active")
        
        // Configuration
        .def("setFragmentationThreshold", &ActivePool::setFragmentationThreshold,
             py::arg("threshold"),
             "Set fragmentation threshold for automatic compaction")
        .def("getFragmentationThreshold", &ActivePool::getFragmentationThreshold,
             "Get fragmentation threshold")
        
        // Capacity
        .def("getMaxAtoms", &ActivePool::getMaxAtoms,
             "Get maximum atom capacity")
        .def("getMaxResidues", &ActivePool::getMaxResidues,
             "Get maximum residue capacity")
        .def("canInsert", &ActivePool::canInsert,
             py::arg("atomCount"),
             "Check if there's space for insertion")
        
        // Statistics
        .def("getStatistics", &ActivePool::getStatistics,
             py::return_value_policy::reference_internal,
             "Get pool statistics")
        .def("resetStatistics", &ActivePool::resetStatistics,
             "Reset statistics")
        
        .def("__repr__", [](const ActivePool& pool) {
            auto counts = pool.getActiveCounts();
            return "<ActivePool atoms=" + std::to_string(counts.first) + "/" + 
                   std::to_string(pool.getMaxAtoms()) + 
                   " residues=" + std::to_string(counts.second) + "/" +
                   std::to_string(pool.getMaxResidues()) + ">";
        });
    
    // Bind ActivePool::Statistics
    py::class_<ActivePool::Statistics>(movement, "ActivePoolStatistics")
        .def_readonly("totalInserts", &ActivePool::Statistics::totalInserts)
        .def_readonly("totalDeletes", &ActivePool::Statistics::totalDeletes)
        .def_readonly("compactions", &ActivePool::Statistics::compactions)
        .def_readonly("peakAtoms", &ActivePool::Statistics::peakAtoms)
        .def_readonly("peakResidues", &ActivePool::Statistics::peakResidues)
        .def_readonly("averageFragmentation", &ActivePool::Statistics::averageFragmentation);
    
    // Add module-level docstring
    movement.doc() = R"pbdoc(
        PyGCMC Movement Module
        ----------------------
        
        This module provides GCMC movement operations including:
        - Insertion/Deletion with cavity bias
        - Translation moves
        - Rotation with configurational bias
        - Active pool memory management
        - Numerical stability through log-space calculations
        - Fragment reservoir management
        
        Example usage:
            import pygcmc
            
            # Create movement module with parameters
            params = pygcmc.movement.MovementParams(temperature=298.15)
            params.chemicalPotential = -15.7  # kJ/mol for water
            params.useCavityBias = True
            params.useConfigBias = True
            
            movement = pygcmc.movement.MovementModule(params)
            
            # Perform movements
            result = movement.attemptInsertion(state)
            if result.accepted:
                print(f"Insertion accepted with ΔE = {result.energyChange} kJ/mol")
    )pbdoc";
    
    // Initialize FragmentReservoir bindings
    bindings::movement::init_fragment_reservoir_bindings(movement);
}

} // namespace simulation
} // namespace bindings
} // namespace pygcmc