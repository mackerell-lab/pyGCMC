// CBMC insertion helper functions for Insertion.cpp
#include "../../energy/EnergyModule.hpp"
// This file contains the CBMC-specific logic that will be integrated into performCavityBiasInsertion

#include "Insertion.hpp"
#include "../pool/ActivePool.hpp"
#include "../common/MovementUtils.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include <cmath>
#include <limits>
#include <random>
#include <chrono>
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;
using namespace utils;

// Helper function to generate trial configurations
std::vector<std::vector<MCAtom>> InsertionMove::generateTrialConfigurations(
    int moleculeType,
    const Vector3& position,
    const MovementParams& params,
    const MCState& state) 
{
    std::vector<std::vector<MCAtom>> trials;
    trials.reserve(params.numConfigTrials);
    
    for (int i = 0; i < params.numConfigTrials; ++i) {
        std::vector<MCAtom> atoms = createMolecule(moleculeType, position);
        
        // Apply random rotation for trials beyond the first
        if (i > 0) {
            // Generate random quaternion
            Quaternion quat = RotationUtils::generateRandomQuaternion();
            double matrix[3][3];
            RotationUtils::quaternionToMatrix(quat, matrix);
            
            // Calculate molecule center
            Vector3 center(0, 0, 0);
            for (const auto& atom : atoms) {
                center.x += atom.x;
                center.y += atom.y;
                center.z += atom.z;
            }
            center.x /= atoms.size();
            center.y /= atoms.size();
            center.z /= atoms.size();
            
            // Rotate atoms around center
            for (auto& atom : atoms) {
                Vector3 rel(atom.x - center.x, atom.y - center.y, atom.z - center.z);
                Vector3 rotated = RotationUtils::rotateVector(rel, matrix);
                atom.x = center.x + rotated.x;
                atom.y = center.y + rotated.y;
                atom.z = center.z + rotated.z;
            }
            
            // Optional: add small translation with PBC
            if (params.configTranslationRange > 0) {
                Vector3 delta = RandomUtils::randomVector(params.configTranslationRange);
                
                // Apply translation and PBC using real box dimensions
                for (auto& atom : atoms) {
                    atom.x += delta.x;
                    atom.y += delta.y;
                    atom.z += delta.z;
                    
                    // Apply PBC with real box from state
                    utils::PBCUtils::applyPBC(atom.x, atom.y, atom.z, 
                                            const_cast<float*>(state.info.box));
                }
            }
        }
        
        trials.push_back(std::move(atoms));
    }
    
    return trials;
}

// Helper function to evaluate trial energies
std::pair<std::vector<double>, int> InsertionMove::evaluateTrialEnergies(
    MCState& state,
    const std::vector<std::vector<MCAtom>>& trials,
    int moleculeType,
    double /*energyBefore*/) 
{
    std::vector<double> deltaEnergies(trials.size());
    int Keff = 0;  // Effective number of valid trials
    
    for (size_t i = 0; i < trials.size(); ++i) {
        // Temporarily insert trial configuration
        int tempResIdx = state.addResidue(MCResidue());
        MCResidue& newRes = state.residues[tempResIdx];
        newRes.atomStart = state.activeAtomCount;
        newRes.atomCount = static_cast<int>(trials[i].size());
        newRes.type = moleculeType;
        newRes.active = true;
        
        for (const auto& atom : trials[i]) {
            state.addAtom(atom);
        }
        
        // Calculate energy of the new molecule with existing system
        // This computes interaction energy without double counting
        platform::cpu::computeMovementEnergyCutoff(state);
        
        // Get energy of new residue (interaction with existing atoms only)
        double energyAfter = state.residues[tempResIdx].energy_vdw 
                           + state.residues[tempResIdx].energy_elec;
        
        // Note: This should be the interaction energy of new molecule with existing system
        // No division by 2 needed as we're only computing new residue's energy
        
        deltaEnergies[i] = energyAfter;  // Since energyBefore should be 0 for new residue
        
        // Check validity
        if (std::isfinite(deltaEnergies[i])) {
            Keff++;
        } else {
            deltaEnergies[i] = std::numeric_limits<double>::infinity();
        }
        
        // Immediately rollback
        state.removeResidue(tempResIdx);
        for (int j = 0; j < newRes.atomCount; ++j) {
            state.removeAtom(state.activeAtomCount - 1);
        }
    }
    
    return {deltaEnergies, Keff};
}

// Helper function to select configuration by Boltzmann weights using Gumbel-max trick
int InsertionMove::selectByBoltzmannWeight(
    const std::vector<double>& deltaEnergies,
    double beta,
    double& logWnew) 
{
    // Calculate log weights: log(exp(-beta*deltaE_i))
    std::vector<double> logWeights(deltaEnergies.size());
    for (size_t i = 0; i < deltaEnergies.size(); ++i) {
        logWeights[i] = -beta * deltaEnergies[i];
    }
    
    // Calculate logWnew for acceptance probability
    logWnew = LogSpaceCalculator::logSumExp(logWeights);
    
    // Use Gumbel-max trick for numerically stable selection
    // Sample Gumbel noise and add to log weights
    std::vector<double> gumbelScores(deltaEnergies.size());
    double maxScore = -std::numeric_limits<double>::infinity();
    int selectedIdx = 0;
    
    for (size_t i = 0; i < deltaEnergies.size(); ++i) {
        // Sample Gumbel(0) noise: -log(-log(U)) where U ~ Uniform(0,1)
        double u = RandomUtils::uniform(1e-10, 1.0);  // Avoid log(0)
        double gumbelNoise = -std::log(-std::log(u));
        
        // Score = log_weight + Gumbel noise
        gumbelScores[i] = logWeights[i] + gumbelNoise;
        
        // Track maximum
        if (gumbelScores[i] > maxScore) {
            maxScore = gumbelScores[i];
            selectedIdx = static_cast<int>(i);
        }
    }
    
    return selectedIdx;
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc