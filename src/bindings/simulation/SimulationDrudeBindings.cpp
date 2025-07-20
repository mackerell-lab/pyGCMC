// src/bindings/simulation/SimulationDrudeBindings.cpp

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "../../simulation/simulation.hpp"
#include "../../model/ModelModule.hpp"
#include "../../platform/cpu/energy/drude/DrudeForce.hpp"

namespace py = pybind11;
using namespace pygcmc;
using namespace pygcmc::platform::cpu;

namespace pygcmc {
namespace bindings {
namespace simulation {

// Global DrudeForce object for simplified API
static std::unique_ptr<DrudeForce> globalDrudeForce;

void init_drude_bindings(py::module& m) {
    // DrudeForce class
    py::class_<DrudeForce>(m, "DrudeForce", "Drude oscillator force calculation")
        .def(py::init<>())
        .def("addParticle", &DrudeForce::addParticle,
            py::arg("drudeIndex"),
            py::arg("parentIndex"),
            py::arg("aniso1Index"),
            py::arg("aniso2Index"),
            py::arg("aniso3Index"),
            py::arg("aniso4Index"),
            py::arg("charge"),
            py::arg("polarizability"),
            py::arg("aniso12"),
            py::arg("aniso34"),
            "Add a Drude particle")
        .def("addScreenedPair", &DrudeForce::addScreenedPair,
            py::arg("dipole1"),
            py::arg("dipole2"),
            py::arg("thole"),
            "Add a screened pair interaction")
        .def("setSCFParameters", &DrudeForce::setSCFParameters,
            py::arg("params"),
            "Set SCF parameters")
        .def("calculateEnergySCF", &DrudeForce::calculateEnergySCF,
            py::arg("state"),
            "Calculate energy using SCF method")
        .def("calculateEnergyOPT3", &DrudeForce::calculateEnergyOPT3,
            py::arg("state"),
            "Calculate energy using OPT3 method")
        .def("setUseOPT3", &DrudeForce::setUseOPT3,
            py::arg("use"),
            "Enable/disable OPT3 algorithm")
        .def("getUseOPT3", &DrudeForce::getUseOPT3,
            "Check if OPT3 is enabled")
        .def("getNumParticles", &DrudeForce::getNumParticles,
            "Get number of Drude particles")
        .def("getNumScreenedPairs", &DrudeForce::getNumScreenedPairs,
            "Get number of screened pairs")
        .def("setAlgorithm", &DrudeForce::setAlgorithm,
            py::arg("algorithm"),
            "Set the algorithm for Drude optimization")
        .def("getAlgorithm", &DrudeForce::getAlgorithm,
            "Get current algorithm")
        .def("setOPT3Coefficients", &DrudeForce::setOPT3Coefficients,
            py::arg("c0"), py::arg("c1"), py::arg("c2"), py::arg("c3"),
            "Set OPT3 coefficients")
        .def("getOPT3Coefficients", &DrudeForce::getOPT3Coefficients,
            "Get current OPT3 coefficients")
        .def("setOPT4Coefficients", &DrudeForce::setOPT4Coefficients,
            py::arg("c0"), py::arg("c1"), py::arg("c2"), py::arg("c3"), py::arg("c4"),
            "Set OPT4 coefficients")
        .def("getOPT4Coefficients", &DrudeForce::getOPT4Coefficients,
            "Get current OPT4 coefficients")
        .def("collectTrainingData", &DrudeForce::collectTrainingData,
            py::arg("state"),
            "Collect training data for OPT3 optimization")
        .def("enableTrainingMode", &DrudeForce::enableTrainingMode,
            py::arg("enable"),
            "Enable/disable training mode")
        .def("isTrainingMode", &DrudeForce::isTrainingMode,
            "Check if training mode is enabled");
    
    // OPT3Coefficients structure
    py::class_<OPT3Coefficients>(m, "OPT3Coefficients", "OPT3 coefficients for Drude model")
        .def(py::init<>())
        .def_readwrite("c0", &OPT3Coefficients::c0, "Zero-order weight")
        .def_readwrite("c1", &OPT3Coefficients::c1, "First-order weight")
        .def_readwrite("c2", &OPT3Coefficients::c2, "Second-order weight")
        .def_readwrite("c3", &OPT3Coefficients::c3, "Third-order weight");
    
    // OPT4Coefficients structure
    py::class_<OPT4Coefficients>(m, "OPT4Coefficients", "OPT4 coefficients for Drude model")
        .def(py::init<>())
        .def_readwrite("c0", &OPT4Coefficients::c0, "Zero-order weight")
        .def_readwrite("c1", &OPT4Coefficients::c1, "First-order weight")
        .def_readwrite("c2", &OPT4Coefficients::c2, "Second-order weight")
        .def_readwrite("c3", &OPT4Coefficients::c3, "Third-order weight")
        .def_readwrite("c4", &OPT4Coefficients::c4, "Fourth-order weight");
    
    // DrudeAlgorithm enum
    py::enum_<DrudeAlgorithm>(m, "DrudeAlgorithm", "Algorithm selection for Drude optimization")
        .value("SCF", DrudeAlgorithm::SCF)
        .value("OPT3", DrudeAlgorithm::OPT3)
        .value("OPT4", DrudeAlgorithm::OPT4)
        .value("AdaptiveOPT", DrudeAlgorithm::AdaptiveOPT)
        .value("HybridOPT", DrudeAlgorithm::HybridOPT)
        .value("SmartOPT3", DrudeAlgorithm::SmartOPT3)
        .value("FBP", DrudeAlgorithm::FBP)
        .value("ConjugateGradient", DrudeAlgorithm::ConjugateGradient);
    
    // Vec3 binding for force vectors
    py::class_<DrudeForce::Vec3>(m, "Vec3", "3D vector")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def_readwrite("x", &DrudeForce::Vec3::x)
        .def_readwrite("y", &DrudeForce::Vec3::y)
        .def_readwrite("z", &DrudeForce::Vec3::z)
        .def("norm", &DrudeForce::Vec3::norm)
        .def("__repr__", [](const DrudeForce::Vec3& v) {
            return "Vec3(" + std::to_string(v.x) + ", " + 
                   std::to_string(v.y) + ", " + 
                   std::to_string(v.z) + ")";
        });
    
    // OPT3TrainingData structure
    py::class_<DrudeForce::OPT3TrainingData>(m, "OPT3TrainingData", "Training data for OPT3 optimization")
        .def(py::init<>())
        .def_readonly("r0", &DrudeForce::OPT3TrainingData::r0, "Zero-order displacements")
        .def_readonly("r1", &DrudeForce::OPT3TrainingData::r1, "First-order displacements")
        .def_readonly("r2", &DrudeForce::OPT3TrainingData::r2, "Second-order displacements")
        .def_readonly("r3", &DrudeForce::OPT3TrainingData::r3, "Third-order displacements")
        .def_readonly("r_scf", &DrudeForce::OPT3TrainingData::r_scf, "SCF converged positions")
        .def_readonly("parentPos", &DrudeForce::OPT3TrainingData::parentPos, "Parent positions");
    
    // DrudeSCFParams structure
    py::class_<DrudeSCFParams>(m, "DrudeSCFParams", "SCF parameters for Drude optimization")
        .def(py::init<>())
        .def_readwrite("tolerance", &DrudeSCFParams::tolerance, 
                      "Force tolerance for convergence (kJ/mol/nm)")
        .def_readwrite("maxIterations", &DrudeSCFParams::maxIterations, 
                      "Maximum SCF iterations")
        .def_readwrite("dampingFactor", &DrudeSCFParams::dampingFactor, 
                      "Damping factor for large forces")
        .def_readwrite("forceCutoff", &DrudeSCFParams::forceCutoff, 
                      "Force cutoff for damping (in units of tolerance)")
        .def_readwrite("maxDrudeDistance", &DrudeSCFParams::maxDrudeDistance,
                      "Maximum Drude-parent distance (nm) - hard wall constraint");
    
    // DrudeParticle structure (for inspection)
    py::class_<DrudeParticle>(m, "DrudeParticle", "Drude oscillator parameters")
        .def_readonly("drudeIndex", &DrudeParticle::drudeIndex)
        .def_readonly("parentIndex", &DrudeParticle::parentIndex)
        .def_readonly("charge", &DrudeParticle::charge)
        .def_readonly("polarizability", &DrudeParticle::polarizability)
        .def_readonly("kIsotropic", &DrudeParticle::kIsotropic);
    
    // ScreenedPair structure
    py::class_<ScreenedPair>(m, "ScreenedPair", "Screened dipole pair")
        .def_readonly("dipole1", &ScreenedPair::dipole1)
        .def_readonly("dipole2", &ScreenedPair::dipole2)
        .def_readonly("thole", &ScreenedPair::thole);
    
    // Simplified global API functions
    m.def("initializeDrudeForce", 
        []() {
            globalDrudeForce = std::make_unique<DrudeForce>();
        },
        "Initialize the global Drude force object");
    
    m.def("clearDrudeForce",
        []() {
            globalDrudeForce.reset();
        },
        "Clear the global Drude force object");
    
    m.def("addDrudeParticle",
        [](int drudeIndex, int parentIndex, 
           int aniso1Index, int aniso2Index,
           int aniso3Index, int aniso4Index,
           double charge, double polarizability,
           double aniso12, double aniso34) {
            if (!globalDrudeForce) {
                throw std::runtime_error("Drude force not initialized. Call initializeDrudeForce() first.");
            }
            return globalDrudeForce->addParticle(drudeIndex, parentIndex,
                                                aniso1Index, aniso2Index,
                                                aniso3Index, aniso4Index,
                                                charge, polarizability,
                                                aniso12, aniso34);
        },
        py::arg("drudeIndex"),
        py::arg("parentIndex"),
        py::arg("aniso1Index") = -1,
        py::arg("aniso2Index") = -1,
        py::arg("aniso3Index") = -1,
        py::arg("aniso4Index") = -1,
        py::arg("charge"),
        py::arg("polarizability"),
        py::arg("aniso12") = 1.0,
        py::arg("aniso34") = 1.0,
        "Add a Drude particle to the force");
    
    m.def("addDrudeScreenedPair",
        [](int dipole1, int dipole2, double thole) {
            if (!globalDrudeForce) {
                throw std::runtime_error("Drude force not initialized. Call initializeDrudeForce() first.");
            }
            globalDrudeForce->addScreenedPair(dipole1, dipole2, thole);
        },
        py::arg("dipole1"),
        py::arg("dipole2"),
        py::arg("thole"),
        "Add a screened pair interaction between two dipoles");
    
    m.def("setDrudeSCFParameters",
        [](const DrudeSCFParams& params) {
            if (!globalDrudeForce) {
                throw std::runtime_error("Drude force not initialized. Call initializeDrudeForce() first.");
            }
            globalDrudeForce->setSCFParameters(params);
        },
        py::arg("params"),
        "Set SCF parameters for Drude optimization");
    
    m.def("setDrudeSCFTolerance",
        [](double tolerance) {
            if (!globalDrudeForce) {
                throw std::runtime_error("Drude force not initialized. Call initializeDrudeForce() first.");
            }
            DrudeSCFParams params;
            params.tolerance = tolerance;
            params.maxIterations = 50;
            params.dampingFactor = 0.5;
            params.forceCutoff = 10.0;
            params.maxDrudeDistance = 0.02;  // Default hard wall
            globalDrudeForce->setSCFParameters(params);
        },
        py::arg("tolerance"),
        "Set SCF convergence tolerance (convenience function)");
    
    m.def("computeSystemEnergyDrude",
        [](model::MCState& state) -> std::tuple<double, py::dict> {
            if (!globalDrudeForce) {
                throw std::runtime_error("Drude force not initialized. Call initializeDrudeForce() first.");
            }
            
            double energy = globalDrudeForce->calculateEnergySCF(state);
            
            // Create energy dictionary
            py::dict energy_dict;
            energy_dict["drude_harmonic"] = energy;  // For now, just the total
            energy_dict["drude_screened"] = 0.0;     // Will be separated later
            energy_dict["total"] = energy;
            
            return std::make_tuple(energy, energy_dict);
        },
        py::arg("state"),
        "Calculate Drude energy using SCF method. Returns (total_energy, energy_dict)");
    
    m.def("getNumDrudeParticles",
        []() {
            if (!globalDrudeForce) {
                return 0;
            }
            return globalDrudeForce->getNumParticles();
        },
        "Get the number of Drude particles");
    
    m.def("getNumDrudeScreenedPairs",
        []() {
            if (!globalDrudeForce) {
                return 0;
            }
            return globalDrudeForce->getNumScreenedPairs();
        },
        "Get the number of screened pairs");
    
    // Add OPT3 control functions
    m.def("setDrudeUseOPT3",
        [](bool useOPT3) {
            if (!globalDrudeForce) {
                throw std::runtime_error("Drude force not initialized. Call initializeDrudeForce() first.");
            }
            globalDrudeForce->setUseOPT3(useOPT3);
        },
        py::arg("useOPT3"),
        "Enable/disable OPT3 algorithm for Drude SCF");
    
    m.def("getDrudeUseOPT3",
        []() {
            if (!globalDrudeForce) {
                return false;
            }
            return globalDrudeForce->getUseOPT3();
        },
        "Check if OPT3 is enabled for Drude SCF");
    
    // Create Drude SWM4-NDP water model helper
    m.def("createDrudeSWM4Water",
        [](model::MCState& state, int startIndex) {
            if (!globalDrudeForce) {
                throw std::runtime_error("Drude force not initialized. Call initializeDrudeForce() first.");
            }
            
            // Add Drude particle for oxygen
            // SWM4-NDP parameters
            double charge = -1.71636;
            const double ONE_4PI_EPS0 = 138.935456;  // kJ/mol·nm·e^-2
            double polarizability = ONE_4PI_EPS0 * 1.71636 * 1.71636 / (100000 * 4.184);
            
            return globalDrudeForce->addParticle(
                startIndex + 1,  // Drude index
                startIndex,      // Parent (oxygen) index
                -1, -1, -1, -1,  // No anisotropy
                charge,
                polarizability,
                1.0, 1.0         // Isotropic
            );
        },
        py::arg("state"),
        py::arg("startIndex"),
        "Add Drude particle for SWM4-NDP water model. Assumes atoms are ordered: O, D, H1, H2, M");
}

void init_drude_bindings(py::module& m);

} // namespace simulation
} // namespace bindings
} // namespace pygcmc