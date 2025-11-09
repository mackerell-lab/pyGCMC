#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_RESULT_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_RESULT_HPP

#include <string>
#include "GrandCanonicalTerms.hpp"

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
    double rosenbluthWeight = 1.0;        // Rosenbluth weight (W/K) for CBMC diagnostics
    int numConfigTrials = 1;              // Number of configuration trials used
    int cbmcTrialsUsed = 1;              // Effective CBMC trials contributing to acceptance
    
    // Proposal layer fields (filled when USE_PROPOSAL_LAYER is enabled)
    bool usedCavity = false;              // Whether cavity was used in proposal
    int mproposal = -1;                   // Number of proposal positions (-1 if not applicable)
    double vregion = -1.0;                 // Region volume in nm³ (-1 if not applicable)
    
    // Enhanced diagnostics (P2 strengthening)
    int proposalMode = -1;                // Current proposal mode (0-4, -1 if not set)
    int selectedType = -1;                // ProposalType enum value actually used
    double proposalTimeMs = -1.0;         // Time for proposal generation (ms)
    double findCavTimeMs = -1.0;          // Time for cavity finding (ms)
    
    // Latest proposal info (optional, filled when fillProposalInfo=true)
    double proposalNorm = -1.0;           // q_norm for detailed balance diagnostics
    double proposalPosX = 0.0;            // Proposed position X (nm) - same unit as ProposalInfo
    double proposalPosY = 0.0;            // Proposed position Y (nm) - same unit as ProposalInfo
    double proposalPosZ = 0.0;            // Proposed position Z (nm) - same unit as ProposalInfo
    bool proposalInfoFilled = false;      // Whether proposal info was filled

    // Detailed acceptance metrics
    double logAcceptanceRatio = 0.0;       // log(r) before clamp
    double logProposalForward = 0.0;       // log(p_forward) from scheduler layer
    double logProposalReverse = 0.0;       // log(p_reverse) from scheduler layer
    double volumeNm3 = 0.0;               // Total box volume used in acceptance
    double effectiveVolumeNm3 = 0.0;      // Effective volume after cavity weighting
    double cavityVolumeNm3 = 0.0;         // Actual cavity volume (if available)
    double logVolume = 0.0;               // log(volumeNm3)
    double logCavityFactor = 0.0;         // log(cavityBiasFactor)
    double lambdaNm = 1.0;                // Thermal de Broglie wavelength
    double logLambda3 = 0.0;              // 3*log(lambdaNm)
    double logWForward = 0.0;             // log Rosenbluth forward weight (insertion)
    double logWReverse = 0.0;             // log Rosenbluth reverse weight (deletion)
    
    // Performance metrics
    double computeTimeMs = 0.0;           // Time taken for the move in milliseconds
    
    // Optional detailed information
    std::string rejectReason;             // Reason for rejection (if applicable)
    bool numericalError = false;          // Flag for numerical issues

    // Unified grand-canonical bookkeeping
    bool hasGrandTerms = false;                          // Whether grandTerms/logs are valid
    gcmc::GrandCanonicalTerms grandTerms;                // Logged factors for this move
    gcmc::GrandCanonicalEvaluation grandEvaluation;      // Cached evaluation result
    
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
