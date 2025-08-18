#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "platform/cpu/movement/MovementModule.hpp"
#include "model/montecarlo/MCMain.hpp"  // For MCState
// All necessary types are included via MovementModule.hpp

namespace py = pybind11;
using namespace pygcmc::platform::cpu::movement;

namespace pygcmc {
namespace bindings {
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
        .def("updateDerivedParameters", &MovementParams::updateDerivedParameters,
                      "Update derived parameters after changing temperature")
        .def("__repr__", [](const MovementParams& p) {
            return "<MovementParams T=" + std::to_string(p.temperature) + 
                   " μ=" + std::to_string(p.chemicalPotential) + ">";
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
}

} // namespace simulation
} // namespace bindings
} // namespace pygcmc