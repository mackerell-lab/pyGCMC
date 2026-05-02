#include "DrudeSequentialOptimizer.hpp"
#include "platform/cpu/energy/drude/DrudeFastFBP.hpp"
#include "platform/cpu/energy/drude/DrudeMultiStage.hpp"
#include <sstream>
#include <algorithm>
#include <stdexcept>

namespace drude {

using DrudeAlgorithm = pygcmc::platform::cpu::DrudeAlgorithm;
using MCState = pygcmc::model::montecarlo::MCState;

DrudeSequentialOptimizer::DrudeSequentialOptimizer(const OptimizationSequence& sequence)
    : sequence_(sequence) {}

std::unique_ptr<DrudeSequentialOptimizer> DrudeSequentialOptimizer::fromPythonKwargs(
    const std::vector<std::pair<std::string, double>>& kwargs) {
    auto optimizer = std::make_unique<DrudeSequentialOptimizer>();

    // Add steps in the order they appear in kwargs
    for (const auto& [algorithm, parameter] : kwargs) {
        // Skip if parameter is 0 or negative (treated as disabled)
        if (parameter > 0) {
            optimizer->addStep(algorithm, parameter);
        }
    }

    return optimizer;
}

double DrudeSequentialOptimizer::optimize(MCState& state) {
    if (sequence_.empty()) {
        throw std::runtime_error("DrudeSequentialOptimizer: No optimization steps defined");
    }

    double energy = 0.0;

    // Execute each step in sequence
    for (size_t i = 0; i < sequence_.size(); ++i) {
        const auto& [algorithm_name, parameter] = sequence_[i];

        try {
            DrudeAlgorithm algorithm = stringToAlgorithm(algorithm_name);
            energy = executeStep(state, algorithm, parameter);

            // Log step completion (could be made optional)
            // std::cout << "Step " << (i+1) << "/" << sequence_.size()
            //           << ": " << algorithm_name << "(" << parameter << ")"
            //           << " -> E = " << energy << " kJ/mol" << std::endl;

        } catch (const std::exception& e) {
            std::stringstream ss;
            ss << "Error in step " << (i+1) << " (" << algorithm_name << "): " << e.what();
            throw std::runtime_error(ss.str());
        }
    }

    return energy;
}

void DrudeSequentialOptimizer::addStep(const std::string& algorithm, double parameter) {
    sequence_.emplace_back(algorithm, parameter);
}

void DrudeSequentialOptimizer::clearSequence() {
    sequence_.clear();
}

std::string DrudeSequentialOptimizer::toString() const {
    if (sequence_.empty()) {
        return "DrudeOptimizer[empty]";
    }

    std::stringstream ss;
    ss << "DrudeOptimizer[";

    for (size_t i = 0; i < sequence_.size(); ++i) {
        const auto& [algorithm, parameter] = sequence_[i];

        if (i > 0) {
            ss << " → ";
        }

        ss << algorithm << "(";

        // Format parameter based on algorithm type
        if (algorithm == "scf") {
            ss << parameter;  // Tolerance, show as float
        } else {
            ss << static_cast<int>(parameter);  // Iterations, show as int
        }

        ss << ")";
    }

    ss << "]";
    return ss.str();
}

DrudeAlgorithm DrudeSequentialOptimizer::stringToAlgorithm(const std::string& name) const {
    // Convert to lowercase for case-insensitive matching
    std::string lower_name = name;
    std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(), ::tolower);

    // Map string names to algorithm enums
    static const std::unordered_map<std::string, DrudeAlgorithm> algorithmMap = {
        {"direct", DrudeAlgorithm::Direct},
        {"fast_fbp", DrudeAlgorithm::FastFBP},
        {"fastfbp", DrudeAlgorithm::FastFBP},
        {"fbp", DrudeAlgorithm::FastFBP},
        {"tcg", DrudeAlgorithm::TCG},
        {"scf", DrudeAlgorithm::SCF},
        {"hybrid", DrudeAlgorithm::Hybrid},
        {"multistage", DrudeAlgorithm::MultiStage},
        {"multi_stage", DrudeAlgorithm::MultiStage}
    };

    auto it = algorithmMap.find(lower_name);
    if (it != algorithmMap.end()) {
        return it->second;
    }

    throw std::invalid_argument("Unknown algorithm: " + name);
}

double DrudeSequentialOptimizer::executeStep(MCState& state, DrudeAlgorithm algorithm, double parameter) {
    // Configure algorithm-specific parameters
    configureAlgorithm(algorithm, parameter);

    // Execute the algorithm
    return pygcmc::platform::cpu::DrudeComplete::calculateEnergy(state, algorithm);
}

void DrudeSequentialOptimizer::configureAlgorithm(DrudeAlgorithm algorithm, double parameter) {
    using namespace pygcmc::platform::cpu;

    switch (algorithm) {
        case DrudeAlgorithm::SCF: {
            // Configure SCF tolerance
            DrudeSCFParams params;
            params.tolerance = parameter;
            params.maxIterations = 100;  // Default max iterations
            DrudeComplete::setParameters(params);
            break;
        }

        case DrudeAlgorithm::FastFBP: {
            // Configure FastFBP iterations
            auto* config = DrudeComplete::getFastFBPOptimizer();
            if (config) {
                config->setIterations(static_cast<int>(parameter));
            }
            break;
        }

        case DrudeAlgorithm::TCG: {
            // Configure TCG iterations through SCF params
            // TCG uses params.maxIterations for iteration count
            DrudeSCFParams params;
            params.maxIterations = static_cast<int>(parameter);
            // Keep other params at default values
            params.tolerance = 10.0;
            params.dampingFactor = 0.5;
            DrudeComplete::setParameters(params);
            break;
        }

        case DrudeAlgorithm::Direct:
        case DrudeAlgorithm::Hybrid:
        case DrudeAlgorithm::MultiStage:
            // These algorithms don't require parameter configuration
            // or have their own internal configuration
            break;

        default:
            break;
    }
}

// DrudeOptimizerBuilder implementation

DrudeOptimizerBuilder& DrudeOptimizerBuilder::direct(int iterations) {
    return add("direct", static_cast<double>(iterations));
}

DrudeOptimizerBuilder& DrudeOptimizerBuilder::fastFBP(int iterations) {
    return add("fast_fbp", static_cast<double>(iterations));
}

DrudeOptimizerBuilder& DrudeOptimizerBuilder::tcg(int iterations) {
    return add("tcg", static_cast<double>(iterations));
}

DrudeOptimizerBuilder& DrudeOptimizerBuilder::scf(double tolerance) {
    return add("scf", tolerance);
}


DrudeOptimizerBuilder& DrudeOptimizerBuilder::add(const std::string& algorithm, double parameter) {
    sequence_.emplace_back(algorithm, parameter);
    return *this;
}

std::unique_ptr<DrudeSequentialOptimizer> DrudeOptimizerBuilder::build() {
    return std::make_unique<DrudeSequentialOptimizer>(sequence_);
}

} // namespace drude
