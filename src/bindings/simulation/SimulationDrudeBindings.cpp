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
                       "Force tolerance for convergence (kJ/mol/nm)")
        .def_readwrite("maxIterations", &DrudeSCFParams::maxIterations,
                       "Maximum number of SCF iterations")
        .def_readwrite("dampingFactor", &DrudeSCFParams::dampingFactor,
                       "Damping factor for stability (0-1)")
        .def_readwrite("maxDrudeDistance", &DrudeSCFParams::maxDrudeDistance,
                       "Maximum allowed Drude-parent distance (nm)")
        .def_readwrite("enableHardWall", &DrudeSCFParams::enableHardWall,
                       "Enable hard wall constraint (default: False, matching OpenMM/CHARMM)");
    
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
        .value("Hybrid", DrudeAlgorithm::Hybrid, "Hybrid FastFBP+SCF strategy");
    
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
}

} // namespace simulation
} // namespace bindings
} // namespace pygcmc