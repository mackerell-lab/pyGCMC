#ifndef DRUDE_SEQUENTIAL_OPTIMIZER_HPP
#define DRUDE_SEQUENTIAL_OPTIMIZER_HPP

#include <vector>
#include <string>
#include <utility>
#include <memory>
#include <unordered_map>
#include <functional>
#include "platform/cpu/energy/drude/DrudeStructures.hpp"
#include "platform/cpu/energy/drude/DrudeMain.hpp"

// Forward declaration of MCState
namespace pygcmc {
namespace model {
namespace montecarlo {
    class MCState;
}
}
}

namespace drude {

// Type aliases
using DrudeAlgorithm = pygcmc::platform::cpu::DrudeAlgorithm;
using MCState = pygcmc::model::montecarlo::MCState;

/**
 * DrudeSequentialOptimizer - Order-sensitive Drude optimizer
 *
 * This optimizer executes Drude optimization algorithms in the exact order specified,
 * allowing DrudeOptimizer(direct=1, fast_fbp=5) to behave differently from
 * DrudeOptimizer(fast_fbp=5, direct=1).
 *
 * Example usage from Python:
 *   opt = DrudeOptimizer(direct=1, fast_fbp=5, tcg=3)
 *   energy = opt.optimize(state)
 *
 * The optimization sequence will be: Direct(1) → FastFBP(5) → TCG(3)
 */
class DrudeSequentialOptimizer {
public:
    using AlgorithmStep = std::pair<std::string, double>;
    using OptimizationSequence = std::vector<AlgorithmStep>;

    /**
     * Default constructor - creates empty sequence
     */
    DrudeSequentialOptimizer() = default;

    /**
     * Construct from ordered sequence of algorithm steps
     * @param sequence Vector of (algorithm_name, parameter) pairs
     */
    explicit DrudeSequentialOptimizer(const OptimizationSequence& sequence);

    /**
     * Factory method for Python bindings that preserves kwargs order
     * @param kwargs Python keyword arguments in order
     * @return Unique pointer to new optimizer instance
     */
    static std::unique_ptr<DrudeSequentialOptimizer> fromPythonKwargs(
        const std::vector<std::pair<std::string, double>>& kwargs);

    /**
     * Execute the optimization sequence on the given state
     * @param state The molecular state to optimize
     * @return Final energy after optimization
     */
    double optimize(MCState& state);

    /**
     * Add an algorithm step to the sequence
     * @param algorithm Algorithm name (direct, fast_fbp, tcg, scf, etc.)
     * @param parameter Algorithm-specific parameter (iterations or tolerance)
     */
    void addStep(const std::string& algorithm, double parameter);

    /**
     * Clear the optimization sequence
     */
    void clearSequence();

    /**
     * Get the current optimization sequence
     * @return Vector of algorithm steps
     */
    const OptimizationSequence& getSequence() const { return sequence_; }

    /**
     * Get a string representation of the optimization sequence
     * @return String like "Direct(1) → FastFBP(5) → TCG(3)"
     */
    std::string toString() const;

private:
    OptimizationSequence sequence_;

    /**
     * Convert algorithm name to DrudeAlgorithm enum
     * @param name Algorithm name string
     * @return Corresponding DrudeAlgorithm enum value
     */
    DrudeAlgorithm stringToAlgorithm(const std::string& name) const;

    /**
     * Execute a single algorithm step
     * @param state The molecular state
     * @param algorithm The algorithm to execute
     * @param parameter Algorithm-specific parameter
     * @return Energy after this step
     */
    double executeStep(MCState& state, DrudeAlgorithm algorithm, double parameter);

    /**
     * Configure algorithm-specific parameters before execution
     * @param algorithm The algorithm to configure
     * @param parameter The parameter value
     */
    void configureAlgorithm(DrudeAlgorithm algorithm, double parameter);
};

/**
 * Builder pattern for creating DrudeSequentialOptimizer
 * Provides a fluent API for constructing optimization sequences
 */
class DrudeOptimizerBuilder {
public:
    DrudeOptimizerBuilder() = default;

    /**
     * Add Direct polarization step
     * @param iterations Number of iterations (usually 1)
     * @return Reference to this builder
     */
    DrudeOptimizerBuilder& direct(int iterations = 1);

    /**
     * Add Fast Force Balance Predictor step
     * @param iterations Number of iterations (default 5)
     * @return Reference to this builder
     */
    DrudeOptimizerBuilder& fastFBP(int iterations = 5);

    /**
     * Add Truncated Conjugate Gradient step
     * @param iterations Number of iterations (default 3)
     * @return Reference to this builder
     */
    DrudeOptimizerBuilder& tcg(int iterations = 3);

    /**
     * Add Self-Consistent Field step
     * @param tolerance Convergence tolerance in nm (default 0.01)
     * @return Reference to this builder
     */
    DrudeOptimizerBuilder& scf(double tolerance = 0.01);


    /**
     * Add custom algorithm step
     * @param algorithm Algorithm name
     * @param parameter Algorithm parameter
     * @return Reference to this builder
     */
    DrudeOptimizerBuilder& add(const std::string& algorithm, double parameter);

    /**
     * Build the final optimizer
     * @return Unique pointer to the constructed optimizer
     */
    std::unique_ptr<DrudeSequentialOptimizer> build();

private:
    DrudeSequentialOptimizer::OptimizationSequence sequence_;
};

} // namespace drude

#endif // DRUDE_SEQUENTIAL_OPTIMIZER_HPP
