/**
 * @file SimulationDrudeBindings.cpp
 * @brief Python bindings for Drude oscillator module
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "platform/cpu/energy/drude/DrudeMain.hpp"
#include "platform/cpu/energy/drude/DrudeStructures.hpp"
#include "platform/cpu/energy/drude/DrudeFastFBP.hpp"
#include "platform/cpu/energy/drude/DrudeHybrid.hpp"
#include "platform/cpu/energy/drude/DrudeMultiStage.hpp"
#include "platform/cpu/energy/drude/DrudeSequentialOptimizer.hpp"
#include "platform/cpu/energy/drude/exp/DrudeExperimentalCore.hpp"
#include "platform/cpu/energy/drude/exp/DrudeSCFOM.hpp"

namespace py = pybind11;

namespace pygcmc {
namespace bindings {
namespace simulation {

void init_drude_bindings(py::module& m) {
    using namespace pygcmc::platform::cpu;
    // DrudeConstants
    py::module constants = m.def_submodule("DrudeConstants", "Physical constants for Drude calculations");
    constants.attr("ONE_4PI_EPS0") = DrudeConstants::ONE_4PI_EPS0;
    constants.attr("DRUDE_MASS") = DrudeConstants::DRUDE_MASS;
    
    // DrudeParticle
    py::class_<DrudeParticle>(m, "DrudeParticle", "Drude oscillator parameters for a single particle")
        .def(py::init<>())
        .def_readwrite("drudeIndex", &DrudeParticle::drudeIndex)
        .def_readwrite("parentIndex", &DrudeParticle::parentIndex)
        .def_readwrite("aniso1Index", &DrudeParticle::aniso1Index)
        .def_readwrite("aniso2Index", &DrudeParticle::aniso2Index)
        .def_readwrite("aniso3Index", &DrudeParticle::aniso3Index)
        .def_readwrite("aniso4Index", &DrudeParticle::aniso4Index)
        .def_readwrite("charge", &DrudeParticle::charge)
        .def_readwrite("polarizability", &DrudeParticle::polarizability)
        .def_readwrite("aniso12", &DrudeParticle::aniso12)
        .def_readwrite("aniso34", &DrudeParticle::aniso34)
        .def_readonly("kSpring", &DrudeParticle::kSpring)
        .def_readonly("kAniso1", &DrudeParticle::kAniso1)
        .def_readonly("kAniso2", &DrudeParticle::kAniso2)
        .def("computeSpringConstants", &DrudeParticle::computeSpringConstants,
             "Compute derived spring constants from charge and polarizability");
    
    // ScreenedPair
    py::class_<ScreenedPair>(m, "ScreenedPair", "Thole-screened dipole-dipole interaction")
        .def(py::init<>())
        .def_readwrite("dipole1", &ScreenedPair::dipole1)
        .def_readwrite("dipole2", &ScreenedPair::dipole2)
        .def_readwrite("thole", &ScreenedPair::thole);
    
    // DrudeSCFParams
    py::class_<DrudeSCFParams>(m, "DrudeSCFParams", "SCF convergence parameters")
        .def(py::init<>())
        .def_readwrite("tolerance", &DrudeSCFParams::tolerance,
                       "Residual tolerance for convergence (kJ/mol/nm)")
        .def_readwrite("displacementTolerance", &DrudeSCFParams::displacementTolerance,
                       "Displacement tolerance for convergence (nm)")
        .def_readwrite("maxIterations", &DrudeSCFParams::maxIterations,
                       "Maximum number of SCF iterations")
        .def_readwrite("dampingFactor", &DrudeSCFParams::dampingFactor,
                       "Initial damping factor for stability (0-1)")
        .def_readwrite("maxStep", &DrudeSCFParams::maxStep,
                       "Maximum step size per iteration (nm)")
        .def_readwrite("maxDrudeDistance", &DrudeSCFParams::maxDrudeDistance,
                       "Maximum allowed Drude-parent distance (nm)")
        .def_readwrite("diisStartIter", &DrudeSCFParams::diisStartIter,
                       "Iteration to start DIIS acceleration")
        .def_readwrite("diisMaxHistory", &DrudeSCFParams::diisMaxHistory,
                       "Maximum DIIS history size")
        .def_readwrite("enableHardWall", &DrudeSCFParams::enableHardWall,
                       "Enable hard wall constraint (default: False, matching OpenMM/CHARMM)")
        .def_readwrite("excludePartnerParentInExternalField", &DrudeSCFParams::excludePartnerParentInExternalField,
                       "Exclude partner partners in external field to avoid double counting (default: False)")
        .def_readwrite("includeCoulombEnergy", &DrudeSCFParams::includeCoulombEnergy,
                       "Include Coulomb energy in calculateEnergy (default: False, spring-only for production)")
        .def_readwrite("enableAdaptiveDamping", &DrudeSCFParams::enableAdaptiveDamping,
                       "Enable adaptive damping based on spectral radius (default: False)")
        .def_readwrite("requireConvergence", &DrudeSCFParams::requireConvergence,
                       "Throw exception if SCF does not converge (default: False)")
        .def_readwrite("logLevel", &DrudeSCFParams::logLevel,
                       "Logging level: 0=silent, 1=brief, 2=verbose")
        .def_readwrite("tholeMode", &DrudeSCFParams::tholeMode,
                       "Thole screening mode: StandardS1 or OpenMMCompat (default: StandardS1)")
        .def_readwrite("compatSmallUSoftening", &DrudeSCFParams::compatSmallUSoftening,
                       "Enable small-u softening for P-D interactions (OpenMMCompat only)")
        .def_readwrite("compatUSoftenStart", &DrudeSCFParams::compatUSoftenStart,
                       "u0: below this, S1 approaches 1 (default: 0.2)")
        .def_readwrite("compatUSoftenEnd", &DrudeSCFParams::compatUSoftenEnd,
                       "u1: above this, use standard S1 (default: 0.9)");
    
    // TholeMode enum
    py::enum_<TholeMode>(m, "TholeMode", "Thole screening modes")
        .value("StandardS1", TholeMode::StandardS1, "Standard theoretical S1 function")
        .value("OpenMMCompat", TholeMode::OpenMMCompat, "OpenMM-compatible implementation");
    
    // DrudeAlgorithm enum
    py::enum_<DrudeAlgorithm>(m, "DrudeAlgorithm", "Available Drude optimization algorithms")
        .value("SCF", DrudeAlgorithm::SCF, "Self-Consistent Field iteration")
        .value("OPT3", DrudeAlgorithm::OPT3, "3rd order perturbation theory")
        .value("FBP", DrudeAlgorithm::FBP, "Force Balance Predictor")
        .value("FastFBP", DrudeAlgorithm::FastFBP, "Fast Force Balance Predictor (5% accuracy for GCMC)")
        .value("TCG", DrudeAlgorithm::TCG, "Truncated Conjugate Gradient (fixed iterations)")
        .value("TCGv2", DrudeAlgorithm::TCGv2, "Improved TCG with preconditioning")
        .value("LBFGS", DrudeAlgorithm::LBFGS, "L-BFGS optimization (matches OpenMM precision)")
        .value("Direct", DrudeAlgorithm::Direct, "Direct polarization (ignores induced-induced)")
        .value("Hybrid", DrudeAlgorithm::Hybrid, "Hybrid FastFBP+SCF strategy")
        .value("MultiStage", DrudeAlgorithm::MultiStage, "Multi-stage Direct→FastFBP→TCG→SCF optimization");
    
    // OPT3Coefficients
    py::class_<OPT3Coefficients>(m, "OPT3Coefficients", "OPT3 expansion coefficients")
        .def(py::init<>())
        .def_readwrite("c0", &OPT3Coefficients::c0, "Zero-order coefficient")
        .def_readwrite("c1", &OPT3Coefficients::c1, "First-order coefficient")
        .def_readwrite("c2", &OPT3Coefficients::c2, "Second-order coefficient")
        .def_readwrite("c3", &OPT3Coefficients::c3, "Third-order coefficient");
    
    // DrudeComplete static interface
    py::class_<DrudeComplete>(m, "DrudeComplete", "Main interface for Drude force calculations")
        .def_static("calculateEnergy", 
                    py::overload_cast<pygcmc::model::MCState&>(&DrudeComplete::calculateEnergy),
                    py::arg("state"),
                    "Calculate Drude energy with SCF optimization")
        .def_static("calculateEnergy",
                    py::overload_cast<pygcmc::model::MCState&, DrudeAlgorithm>(&DrudeComplete::calculateEnergy),
                    py::arg("state"), py::arg("algorithm"),
                    "Calculate Drude energy using specific algorithm")
        .def_static("setParameters", &DrudeComplete::setParameters,
                    py::arg("params"),
                    "Set global Drude SCF parameters")
        .def_static("addParticle", &DrudeComplete::addParticle,
                    py::arg("particle"),
                    "Add a Drude particle to the system")
        .def_static("addScreenedPair", &DrudeComplete::addScreenedPair,
                    py::arg("pair"),
                    "Add a Thole-screened pair interaction")
        .def_static("clear", &DrudeComplete::clear,
                    "Clear all Drude particles and pairs")
        .def_static("getNumParticles", &DrudeComplete::getNumParticles,
                    "Get the number of Drude particles")
        .def_static("enableASPC", &DrudeComplete::enableASPC,
                    py::arg("enable"),
                    "Enable/disable ASPC history prediction for faster SCF convergence")
        .def_static("isASPCEnabled", &DrudeComplete::isASPCEnabled,
                    "Check if ASPC history prediction is enabled")
        .def_static("clearHistory", &DrudeComplete::clearHistory,
                    "Clear ASPC history (useful when system changes significantly)");
    
    // Free function for Thole screening
    m.def("computeTholeScreening", &computeTholeScreening,
          py::arg("r"), py::arg("alpha_i"), py::arg("alpha_j"), py::arg("thole"),
          "Compute Thole screening function value");
    
    // FastFBP IterationMode enum
    py::enum_<DrudeFastFBP::IterationMode>(m, "FastFBPIterationMode", "Iteration modes for FastFBP")
        .value("Fixed", DrudeFastFBP::IterationMode::Fixed, "Use fixed number of iterations")
        .value("Dynamic", DrudeFastFBP::IterationMode::Dynamic, "Dynamic iterations based on convergence")
        .value("Adaptive", DrudeFastFBP::IterationMode::Adaptive, "Adaptive with parameter adjustment");
    
    // FastFBP ConvergenceStats
    py::class_<DrudeFastFBP::ConvergenceStats>(m, "FastFBPConvergenceStats", "Convergence statistics for FastFBP")
        .def_readonly("actualIterations", &DrudeFastFBP::ConvergenceStats::actualIterations)
        .def_readonly("finalError", &DrudeFastFBP::ConvergenceStats::finalError)
        .def_readonly("convergenceRate", &DrudeFastFBP::ConvergenceStats::convergenceRate)
        .def_readonly("converged", &DrudeFastFBP::ConvergenceStats::converged);
    
    // FastFBPConfig class for configuration
    py::class_<DrudeFastFBP>(m, "FastFBPConfig", "Configuration for FastFBP algorithm")
        .def_static("getInstance", []() -> DrudeFastFBP* {
            return DrudeComplete::getFastFBPOptimizer();
        }, py::return_value_policy::reference, "Get FastFBP optimizer instance")
        .def("setIterations", &DrudeFastFBP::setIterations,
             py::arg("iterations"), "Set number of FBP iterations for fixed mode")
        .def("setCutoff", &DrudeFastFBP::setCutoff,
             py::arg("cutoff"), "Set interaction cutoff for Drude-Drude (nm)")
        .def("setIterationMode", &DrudeFastFBP::setIterationMode,
             py::arg("mode"), "Set iteration mode (Fixed, Dynamic, or Adaptive)")
        .def("setConvergenceTolerance", &DrudeFastFBP::setConvergenceTolerance,
             py::arg("tolerance"), "Set convergence tolerance for dynamic mode (nm)")
        .def("setMaxIterations", &DrudeFastFBP::setMaxIterations,
             py::arg("max_iterations"), "Set maximum iterations for dynamic mode")
        .def("setAdaptiveMode", &DrudeFastFBP::setAdaptiveMode,
             py::arg("enable"), "Enable/disable adaptive parameter adjustment")
        .def("getStats", &DrudeFastFBP::getStats,
             py::return_value_policy::reference_internal,
             "Get convergence statistics from last optimization");
    
    // DrudeHybrid.HybridMode enum
    py::enum_<DrudeHybrid::HybridMode>(m, "HybridMode", "Hybrid optimization strategy modes")
        .value("Fixed", DrudeHybrid::HybridMode::Fixed, "Fixed FastFBP iterations before SCF")
        .value("Dynamic", DrudeHybrid::HybridMode::Dynamic, "Dynamic switching based on convergence")
        .value("Adaptive", DrudeHybrid::HybridMode::Adaptive, "Adaptive with learning");
    
    // DrudeHybrid.HybridStats
    py::class_<DrudeHybrid::HybridStats>(m, "HybridStats", "Statistics from hybrid optimization")
        .def_readonly("fbpIterations", &DrudeHybrid::HybridStats::fbpIterations, "Actual FastFBP iterations used")
        .def_readonly("scfIterations", &DrudeHybrid::HybridStats::scfIterations, "SCF iterations after switching")
        .def_readonly("fbpTime", &DrudeHybrid::HybridStats::fbpTime, "Time spent in FastFBP (seconds)")
        .def_readonly("scfTime", &DrudeHybrid::HybridStats::scfTime, "Time spent in SCF (seconds)")
        .def_readonly("switchError", &DrudeHybrid::HybridStats::switchError, "Error when switching to SCF")
        .def_readonly("converged", &DrudeHybrid::HybridStats::converged, "Final convergence status");
    
    // HybridConfig class for configuration
    py::class_<DrudeHybrid>(m, "HybridConfig", "Configuration for Hybrid FastFBP+SCF algorithm")
        .def_static("getInstance", []() -> DrudeHybrid* {
            auto& drudeCore = DrudeComplete::getDrudeCore();
            return drudeCore.getHybridOptimizer();
        }, py::return_value_policy::reference_internal,
           "Get the Hybrid configuration instance")
        .def("setHybridMode", &DrudeHybrid::setHybridMode,
             py::arg("mode"), "Set hybrid strategy mode")
        .def("setFastFBPIterations", &DrudeHybrid::setFastFBPIterations,
             py::arg("iter"), "Set number of FastFBP iterations before switching")
        .def("setSwitchingThreshold", &DrudeHybrid::setSwitchingThreshold,
             py::arg("threshold"), "Set error threshold for switching to SCF")
        .def("setMaxSCFIterations", &DrudeHybrid::setMaxSCFIterations,
             py::arg("iter"), "Set maximum SCF iterations after switching")
        .def("getStats", &DrudeHybrid::getStats,
             py::return_value_policy::reference_internal,
             "Get statistics from last optimization");
    
    // MultiStageConfig class
    py::class_<DrudeMultiStage::MultiStageConfig>(m, "MultiStageConfig", "Configuration for multi-stage optimization")
        .def(py::init<>())
        // FastFBP stage
        .def_readwrite("minFBPIterations", &DrudeMultiStage::MultiStageConfig::minFBPIterations,
                      "Minimum FastFBP iterations before switching")
        .def_readwrite("maxFBPIterations", &DrudeMultiStage::MultiStageConfig::maxFBPIterations,
                      "Maximum FastFBP iterations")
        .def_readwrite("fbpCutoffFactor", &DrudeMultiStage::MultiStageConfig::fbpCutoffFactor,
                      "Fraction of full cutoff for FastFBP")
        // TCG stage
        .def_readwrite("tcgIterations", &DrudeMultiStage::MultiStageConfig::tcgIterations,
                      "Number of TCG iterations (3 recommended)")
        .def_readwrite("enableTCG", &DrudeMultiStage::MultiStageConfig::enableTCG,
                      "Enable TCG stage")
        .def_readwrite("tcgErrorThreshold", &DrudeMultiStage::MultiStageConfig::tcgErrorThreshold,
                      "Error threshold to switch to TCG")
        // Switching criteria
        .def_readwrite("switchToSCFError", &DrudeMultiStage::MultiStageConfig::switchToSCFError,
                      "Error threshold to switch to SCF")
        .def_readwrite("switchToSCFRate", &DrudeMultiStage::MultiStageConfig::switchToSCFRate,
                      "Convergence rate threshold to switch to SCF")
        // SCF stage
        .def_readwrite("scfTolerance", &DrudeMultiStage::MultiStageConfig::scfTolerance,
                      "Final SCF tolerance")
        .def_readwrite("maxSCFIterations", &DrudeMultiStage::MultiStageConfig::maxSCFIterations,
                      "Maximum SCF iterations")
        .def_readwrite("requireSCF", &DrudeMultiStage::MultiStageConfig::requireSCF,
                      "Require SCF stage (False allows stopping after TCG)")
        // Adaptive parameters
        .def_readwrite("adaptiveMode", &DrudeMultiStage::MultiStageConfig::adaptiveMode,
                      "Enable adaptive parameter adjustment")
        .def_readwrite("densityThreshold", &DrudeMultiStage::MultiStageConfig::densityThreshold,
                      "Density threshold for adaptation (g/cm³)")
        .def_readwrite("polarizabilityThreshold", &DrudeMultiStage::MultiStageConfig::polarizabilityThreshold,
                      "Polarizability threshold for adaptation (nm³)");
    
    // MultiStageStats class
    py::class_<DrudeMultiStage::MultiStageStats>(m, "MultiStageStats", "Statistics from multi-stage optimization")
        // Timing
        .def_readonly("directTime", &DrudeMultiStage::MultiStageStats::directTime, "Direct polarization time (s)")
        .def_readonly("fbpTime", &DrudeMultiStage::MultiStageStats::fbpTime, "FastFBP time (s)")
        .def_readonly("tcgTime", &DrudeMultiStage::MultiStageStats::tcgTime, "TCG time (s)")
        .def_readonly("scfTime", &DrudeMultiStage::MultiStageStats::scfTime, "SCF time (s)")
        // Iterations
        .def_readonly("fbpIterations", &DrudeMultiStage::MultiStageStats::fbpIterations, "FastFBP iterations performed")
        .def_readonly("tcgIterations", &DrudeMultiStage::MultiStageStats::tcgIterations, "TCG iterations performed")
        .def_readonly("scfIterations", &DrudeMultiStage::MultiStageStats::scfIterations, "SCF iterations performed")
        // Errors
        .def_readonly("errorAfterDirect", &DrudeMultiStage::MultiStageStats::errorAfterDirect, "Error after Direct (nm)")
        .def_readonly("errorAfterFBP", &DrudeMultiStage::MultiStageStats::errorAfterFBP, "Error after FastFBP (nm)")
        .def_readonly("errorAfterTCG", &DrudeMultiStage::MultiStageStats::errorAfterTCG, "Error after TCG (nm)")
        .def_readonly("finalError", &DrudeMultiStage::MultiStageStats::finalError, "Final error (nm)")
        // System info
        .def_readonly("systemDensity", &DrudeMultiStage::MultiStageStats::systemDensity, "System density (g/cm³)")
        .def_readonly("avgPolarizability", &DrudeMultiStage::MultiStageStats::avgPolarizability, "Average polarizability (nm³)")
        .def_readonly("converged", &DrudeMultiStage::MultiStageStats::converged, "Final convergence status");
    
    // MultiStage class for configuration
    py::class_<DrudeMultiStage>(m, "MultiStageOptimizer", "Multi-stage Direct→FastFBP→TCG→SCF optimizer")
        .def_static("getInstance", []() -> DrudeMultiStage* {
            return DrudeComplete::getMultiStageOptimizer();
        }, py::return_value_policy::reference,
           "Get the MultiStage optimizer instance")
        .def("setConfig", &DrudeMultiStage::setConfig,
             py::arg("config"), "Set configuration parameters")
        .def("getConfig", &DrudeMultiStage::getConfig,
             py::return_value_policy::reference_internal,
             "Get current configuration")
        .def("getStats", &DrudeMultiStage::getStats,
             py::return_value_policy::reference_internal,
             "Get statistics from last optimization")
        .def("enableDirectStage", &DrudeMultiStage::enableDirectStage,
             py::arg("enable"), "Enable/disable Direct polarization stage")
        .def("enableFBPStage", &DrudeMultiStage::enableFBPStage,
             py::arg("enable"), "Enable/disable FastFBP stage")
        .def("enableTCGStage", &DrudeMultiStage::enableTCGStage,
             py::arg("enable"), "Enable/disable TCG stage")
        .def("enableSCFStage", &DrudeMultiStage::enableSCFStage,
             py::arg("enable"), "Enable/disable SCF stage");
    
    // DrudeSequentialOptimizer - Order-sensitive optimizer
    py::class_<drude::DrudeSequentialOptimizer>(m, "DrudeOptimizer")
        .def(py::init<>(), "Create an empty optimizer")
        .def(py::init([](py::kwargs kwargs) {
            // Capture kwargs in order (Python 3.7+ guarantees order)
            std::vector<std::pair<std::string, double>> sequence;
            
            for (auto item : kwargs) {
                std::string key = py::str(item.first);
                double value = 0.0;
                
                // Handle different parameter types
                if (py::isinstance<py::int_>(item.second)) {
                    value = item.second.cast<int>();
                } else if (py::isinstance<py::float_>(item.second)) {
                    value = item.second.cast<double>();
                } else {
                    throw py::type_error("Parameter value must be numeric");
                }
                
                // Only add if value > 0 (0 means disabled)
                if (value > 0) {
                    sequence.push_back({key, value});
                }
            }
            
            return drude::DrudeSequentialOptimizer::fromPythonKwargs(sequence);
        }), R"pbdoc(
            Create a sequential Drude optimizer with ordered algorithm steps.
            
            Parameters are specified as keyword arguments where the order matters:
            - direct: Number of iterations for direct polarization (usually 1)
            - fast_fbp: Number of iterations for Fast Force Balance Predictor
            - tcg: Number of iterations for Truncated Conjugate Gradient
            - scf: Convergence tolerance for Self-Consistent Field (in nm)
            
            Example:
                # Different orders produce different optimization sequences
                opt1 = DrudeOptimizer(direct=1, fast_fbp=5, tcg=3)
                opt2 = DrudeOptimizer(fast_fbp=5, direct=1, tcg=3)
                
                # For GCMC insertion (fast)
                opt_insert = DrudeOptimizer(direct=1, fast_fbp=3)
                
                # For GCMC deletion (accurate)
                opt_delete = DrudeOptimizer(fast_fbp=5, tcg=3)
        )pbdoc")
        
        .def("optimize", &drude::DrudeSequentialOptimizer::optimize,
             py::arg("state"),
             R"pbdoc(
             Execute the optimization sequence on the given state.
             
             Args:
                 state: MCState object to optimize
                 
             Returns:
                 float: Final energy after optimization (kJ/mol)
         )pbdoc")
        
        .def("add_step", &drude::DrudeSequentialOptimizer::addStep,
             py::arg("algorithm"), py::arg("parameter"),
             R"pbdoc(
             Add an algorithm step to the sequence.
             
             Args:
                 algorithm: Algorithm name ('direct', 'fast_fbp', 'tcg', 'scf', etc.)
                 parameter: Algorithm-specific parameter (iterations or tolerance)
         )pbdoc")
        
        .def("clear_sequence", &drude::DrudeSequentialOptimizer::clearSequence,
             "Clear all steps from the optimization sequence")
        
        .def("get_sequence", &drude::DrudeSequentialOptimizer::getSequence,
             "Get the current optimization sequence as a list of (algorithm, parameter) tuples")
        
        .def("__str__", &drude::DrudeSequentialOptimizer::toString)
        
        .def("__repr__", &drude::DrudeSequentialOptimizer::toString)
        
        .def("__len__", [](const drude::DrudeSequentialOptimizer& self) {
            return self.getSequence().size();
        })
        
        .def("__getitem__", [](const drude::DrudeSequentialOptimizer& self, size_t i) {
            const auto& seq = self.getSequence();
            if (i >= seq.size()) {
                throw py::index_error("Index out of range");
            }
            return seq[i];
        });
    
    // DrudeOptimizerBuilder class (optional, for fluent API)
    py::class_<drude::DrudeOptimizerBuilder>(m, "DrudeOptimizerBuilder")
        .def(py::init<>())
        .def("direct", &drude::DrudeOptimizerBuilder::direct,
             py::arg("iterations") = 1,
             "Add Direct polarization step")
        .def("fast_fbp", &drude::DrudeOptimizerBuilder::fastFBP,
             py::arg("iterations") = 5,
             "Add Fast Force Balance Predictor step")
        .def("tcg", &drude::DrudeOptimizerBuilder::tcg,
             py::arg("iterations") = 3,
             "Add Truncated Conjugate Gradient step")
        .def("scf", &drude::DrudeOptimizerBuilder::scf,
             py::arg("tolerance") = 0.01,
             "Add Self-Consistent Field step")
        .def("add", &drude::DrudeOptimizerBuilder::add,
             py::arg("algorithm"), py::arg("parameter"),
             "Add custom algorithm step")
        .def("build", &drude::DrudeOptimizerBuilder::build,
             "Build the final optimizer");
    
    // ============================================================================
    // EXPERIMENTAL DRUDE IMPLEMENTATION (OpenMM-style)
    // ============================================================================
    // WARNING: DrudeExperimental may cause segmentation faults
    // Use DrudeComplete for production. This is for testing/comparison only.
    // Re-enabled temporarily for OpenMM alignment testing
    
    using namespace pygcmc::platform::cpu::exp;
    
    // Static instance for experimental version
    static DrudeExperimentalCore g_experimentalDrude;
    
    py::class_<DrudeExperimentalCore>(m, "DrudeExperimental", 
                                      "Experimental Drude implementation with OpenMM-style algorithm (WARNING: May segfault)")
        .def_static("instance", []() -> DrudeExperimentalCore& {
            return g_experimentalDrude;
        }, py::return_value_policy::reference,
           "Get singleton instance of experimental Drude")
        .def("clear", &DrudeExperimentalCore::clear,
             "Clear all particles and screened pairs")
        .def("addParticle", &DrudeExperimentalCore::addParticle,
             py::arg("particle"),
             "Add a Drude particle")
        .def("addScreenedPair", &DrudeExperimentalCore::addScreenedPair,
             py::arg("pair"),
             "Add a Thole-screened dipole-dipole pair")
        .def("autoScreenPairs", &DrudeExperimentalCore::autoScreenPairs,
             py::arg("thole"), py::arg("cutoff_nm"),
             "Automatically generate all NBTHOLE pairs within cutoff")
        .def("setParameters", &DrudeExperimentalCore::setParameters,
             py::arg("params"),
             "Set SCF parameters")
        .def("calculateEnergy", &DrudeExperimentalCore::calculateEnergy,
             py::arg("state"),
             "Calculate polarization energy (spring energy only)")
        .def("calculateForces", &DrudeExperimentalCore::calculateForces,
             py::arg("state"), py::arg("forces"),
             "Calculate spring forces")
        .def("getNumParticles", &DrudeExperimentalCore::getNumParticles,
             "Get number of Drude particles")
        .def("getNumScreenedPairs", &DrudeExperimentalCore::getNumScreenedPairs,
             "Get number of screened pairs")
        .def("setIncludeCoulomb", &DrudeExperimentalCore::setIncludeCoulomb,
             py::arg("include"),
             "Enable/disable Coulomb energy (for testing only)")
        .def("getIncludeCoulomb", &DrudeExperimentalCore::getIncludeCoulomb,
             "Check if Coulomb energy is included")
        .def("setDrudeAlgorithm", &DrudeExperimentalCore::setDrudeAlgorithm,
             py::arg("algo"),
             "Set algorithm: 0=S1_POINT_CHARGE (CHARMM), 1=S3S5_DIPOLE_TENSOR, 2=S1_DIPOLE_FIELD, 3=DIRECT_COULOMB")
        .def("getDrudeAlgorithm", &DrudeExperimentalCore::getDrudeAlgorithm,
             "Get current algorithm setting")
        .def("getSpectralRadius", &DrudeExperimentalCore::getSpectralRadius,
             py::arg("state"),
             "Get spectral radius of the polarization system for stability analysis")
        .def("getSCFIterationCount", &DrudeExperimentalCore::getSCFIterationCount,
             "Get last SCF iteration count");
}

} // namespace simulation
} // namespace bindings
} // namespace pygcmc