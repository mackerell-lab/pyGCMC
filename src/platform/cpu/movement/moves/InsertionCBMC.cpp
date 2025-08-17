// CBMC insertion helper functions for Insertion.cpp
// This file contains the CBMC-specific logic that will be integrated into performCavityBiasInsertion

#include "Insertion.hpp"
#include "../pool/ActivePool.hpp"
#include "../common/MovementUtils.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include "../../../../simulation/simulation.hpp"
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
    const MovementParams& params) 
{
    static std::mt19937 rng(params.seed != 0 ? params.seed : 
                            std::chrono::steady_clock::now().time_since_epoch().count());
    
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
            
            // Optional: add small translation
            if (params.configTranslationRange > 0) {
                Vector3 delta = RandomUtils::randomVector(params.configTranslationRange);
                for (auto& atom : atoms) {
                    atom.x += delta.x;
                    atom.y += delta.y;
                    atom.z += delta.z;
                }
            }
        }
        
        trials.push_back(atoms);
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
        
        // Calculate energy after insertion
        simulation::Simulation::computeSystemEnergyCutoff(state);
        double energyAfter = 0.0;
        
        // Only calculate energy for the new residue (last one)
        energyAfter = state.residues[tempResIdx].energy_vdw;
        energyAfter += state.residues[tempResIdx].energy_elec;
        
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

// Helper function to select configuration by Boltzmann weights
int InsertionMove::selectByBoltzmannWeight(
    const std::vector<double>& deltaEnergies,
    double beta,
    double& logWnew) 
{
    static std::random_device rd;
    static std::mt19937 rng(rd());
    static std::uniform_real_distribution<> uniform(0.0, 1.0);
    
    // Calculate logWnew = log(sum(exp(-beta*deltaE_i)))
    std::vector<double> logWeights(deltaEnergies.size());
    for (size_t i = 0; i < deltaEnergies.size(); ++i) {
        logWeights[i] = -beta * deltaEnergies[i];
    }
    logWnew = LogSpaceCalculator::logSumExp(logWeights);
    
    // Calculate cumulative probabilities for selection
    std::vector<double> cumProb(deltaEnergies.size());
    double cumSum = 0.0;
    for (size_t i = 0; i < deltaEnergies.size(); ++i) {
        double logP = -beta * deltaEnergies[i] - logWnew;
        cumSum += std::exp(logP);
        cumProb[i] = cumSum;
    }
    
    // Select configuration
    double r = uniform(rng);
    int selectedIdx = static_cast<int>(deltaEnergies.size()) - 1;
    for (size_t i = 0; i < deltaEnergies.size() - 1; ++i) {
        if (r < cumProb[i]) {
            selectedIdx = static_cast<int>(i);
            break;
        }
    }
    
    return selectedIdx;
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc