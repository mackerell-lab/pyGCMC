#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_RESULT_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_RESULT_HPP

#include <string>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

/**
 * Result of a movement attempt
 */
struct MovementResult {
    // Basic result
    bool accepted = false;                 // Whether the move was accepted
    double energyChange = 0.0;            // Change in energy (kJ/mol)
    double acceptanceProbability = 0.0;   // Calculated acceptance probability
    
    // Move details
    std::string moveType;                 // Type of move: "insert", "delete", "translate", "rotate"
    int residueIndex = -1;                // Index of affected residue
    int moleculeType = 0;                 // Type of molecule (for insertion)
    
    // Bias factors (for analysis)
    double cavityBiasFactor = 1.0;        // Cavity bias correction factor
    double configBiasFactor = 1.0;        // Configurational bias correction factor
    int numConfigTrials = 1;              // Number of configuration trials used
    
    // Performance metrics
    double computeTimeMs = 0.0;           // Time taken for the move in milliseconds
    
    // Optional detailed information
    std::string rejectReason;             // Reason for rejection (if applicable)
    bool numericalError = false;          // Flag for numerical issues
    
    // Constructor
    MovementResult() = default;
    
    // Constructor with basic parameters
    MovementResult(bool acc, double dE, double prob, const std::string& type)
        : accepted(acc), energyChange(dE), acceptanceProbability(prob), moveType(type) {}
    
    // Utility function to check if move was successful
    bool isSuccessful() const { return accepted; }
    
    // Get a summary string
    std::string summary() const {
        return moveType + ": " + (accepted ? "accepted" : "rejected") + 
               ", ΔE=" + std::to_string(energyChange) + " kJ/mol" +
               ", P=" + std::to_string(acceptanceProbability);
    }
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_RESULT_HPP