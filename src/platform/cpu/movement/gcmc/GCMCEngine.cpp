#include "GCMCEngine.hpp"
#include "GCMCAcceptance.hpp"
#include "GCMCConfig.hpp"
#include "../../energy/EnergyModule.hpp"
#include "../../energy/common/EnergyDirectCore.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

// Type aliases for clarity
using MCState = model::montecarlo::MCState;
using MCResidue = model::montecarlo::MCResidue;
using MCAtom = model::montecarlo::MCAtom;

// Constructor
GCMCEngine::GCMCEngine()
    : state_(nullptr),
      reservoir_(nullptr),
      cavityManager_(nullptr),
      configBias_(nullptr),
      acceptanceCalculator_(nullptr),
      energyCallback_(std::make_unique<GCMCEnergyCallback>()),
      temperature_(300.0),
      cutoff_(12.0),
      energyMethod_(EnergyMethod::DIRECT),
      rng_(std::random_device{}()),
      uniform_(0.0, 1.0),
      normal_(0.0, 1.0),
      totalMoves_(0),
      acceptedMoves_(0),
      maxTranslationStep_(2.0),
      maxRotationAngleRad_(0.5),
      useCavityBias_(true) {
    energyCache_.valid = false;
    // Configure default energy callback
    energyCallback_->setEnergyMethod(energyMethod_);
    energyCallback_->setParameters(true, true);  // Use cutoff and PBC by default
}

// Destructor
GCMCEngine::~GCMCEngine() {
}

// Initialize
void GCMCEngine::initialize(MCState* state, FragmentReservoir* reservoir) {
    state_ = state;
    reservoir_ = reservoir;
    energyCache_.invalidate();
}

// Set seed - unified for all RNG components
void GCMCEngine::setSeed(unsigned int seed) {
    // Set engine's RNG seed
    rng_.seed(seed);
    lastSeed_ = seed;  // Store for auto-seeding acceptance
    
    // Also set acceptance calculator's RNG seed if present
    if (acceptanceCalculator_) {
        acceptanceCalculator_->setSeed(seed + 1);  // Use different but deterministic seed
    }
    
    // Set reservoir's RNG seed if it has one
    if (reservoir_) {
        // Note: Add setSeed to FragmentReservoir if it needs random operations
        // reservoir_->setSeed(seed + 2);
    }
    
    // Set cavity manager's RNG seed if it has one
    if (cavityManager_) {
        // Note: Add setSeed to CavityManager if it needs random operations
        // cavityManager_->setSeed(seed + 3);
    }
}

// Attempt insertion
GCMCEngine::MoveResult GCMCEngine::attemptInsertion(int typeId) {
    MoveResult result;
    result.type = MoveResult::INSERT;
    result.fragmentType = typeId;
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    // Get template
    FragmentTemplate* tmpl = reservoir_->getTemplate(typeId);
    if (!tmpl) {
        result.accepted = false;
        return result;
    }
    
    // CRITICAL FIX: Get N BEFORE insertion for correct acceptance calculation
    int N_before = reservoir_->getActiveCount(typeId);

    // Determine number of CBMC trials for this fragment type
    int numTrials = 1;
    if (useConfBias_ && typeId < static_cast<int>(cbmcTrialsPerType_.size())) {
        numTrials = cbmcTrialsPerType_[typeId];
    }

    Vector3 position;
    Quaternion orientation;
    double cbmcBias = 1.0;

    if (useConfBias_ && numTrials > 1) {
        // Use CBMC to select configuration
        TrialConfiguration selected = performCBMCInsertion(typeId, numTrials);
        position = selected.position;
        orientation = selected.orientation;
        cbmcBias = selected.weight * numTrials;  // W_new / K
    } else {
        // Original single configuration generation
        position = (useCavityBias_ && cavityManager_) ?
                  generateCavityPosition() : generateRandomPosition();
        applyPeriodicBoundary(position);  // Ensure position is within PBC
        orientation = generateRandomOrientation();
    }

    result.position = position;
    
    // Calculate energy before insertion
    result.energyBefore = calculateSystemEnergy();
    
    // Create instance
    int instanceId = reservoir_->createInstance(typeId, position, orientation);
    if (instanceId < 0) {
        result.accepted = false;
        return result;
    }
    
    // Synchronize MCState with the new instance
    synchronizeStateWithReservoir(instanceId, true);
    
    // Calculate energy change - optimize for DIRECT mode
    if (energyMethod_ == EnergyMethod::DIRECT) {
        // Fast local ΔE calculation - only compute interaction of new residue with system
        cpu::computeResidueEnergyCutoffPBC(*state_, instanceId);
        const auto& residue = state_->residues[instanceId];
        result.deltaE = residue.energy_vdw + residue.energy_elec;
        result.energyBefore = 0.0;  // Not needed for local calculation
        result.energyAfter = result.deltaE;  // For consistency
    } else {
        // Full system energy for Ewald/PME modes
        result.energyAfter = calculateSystemEnergy();
        result.deltaE = result.energyAfter - result.energyBefore;
    }
    
    // Calculate bias (cavity bias * CBMC bias)
    result.bias = calculateInsertionBias(*tmpl, position, orientation) * cbmcBias;
    
    // Calculate acceptance probability using proper GCMC formula
    bool accept = false;
    double prob = 0.0;
    if (acceptanceCalculator_) {
        // CRITICAL FIX: Use N_before for correct detailed balance
        prob = acceptanceCalculator_->calculateInsertionProbability(
            typeId, N_before, result.deltaE, result.bias);
        
        // Only store probability if configured
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        
        accept = acceptanceCalculator_->acceptMove(prob);
    } else {
        // Fallback to simple acceptance (should not be used in production)
        double beta = 1.0 / (8.314e-3 * temperature_);
        double activity = 100.0; // Default activity
        prob = std::min(1.0, (activity * state_->info.volume / (N_before + 1)) * 
                       std::exp(-beta * result.deltaE) * result.bias);
        // Apply same probability storage logic as main branch
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        // Use the calculated prob directly for consistency
        accept = (uniform_(rng_) < prob);
    }
    
    if (accept) {
        result.accepted = true;
        result.residueIndex = instanceId;
        acceptedMoves_++;
        energyCache_.invalidate();
    } else {
        // Remove the instance and revert state
        reservoir_->deleteInstance(instanceId);
        synchronizeStateWithReservoir(instanceId, false);
        result.accepted = false;
    }
    
    totalMoves_++;
    
    // Sample statistics if configured using a lightweight countdown
    if (collectStats_) {
        if (--statsCountdown_ <= 0) {
            int particleCount = reservoir_ ? reservoir_->getActiveCount() : 0;
            double energy = calculateSystemEnergy();
            statistics_.addSample(totalMoves_, particleCount, energy,
                                  getAcceptanceRate(), temperature_, 100.0);
            statsCountdown_ = std::max(1, statsInterval_);
        }
    }
    
    return result;
}

// Attempt deletion
GCMCEngine::MoveResult GCMCEngine::attemptDeletion(int typeId) {
    MoveResult result;
    result.type = MoveResult::DELETE;
    result.fragmentType = typeId;
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    // CRITICAL FIX: Get N BEFORE deletion for correct acceptance calculation
    int N_before = reservoir_->getActiveCount(typeId);
    
    // Check if any instances exist
    if (N_before == 0) {
        result.accepted = false;
        // Set probability semantics: 0.0 when storing, -1.0 when not
        result.acceptanceProbability = shouldStoreProbability() ? 0.0 : -1.0;
        return result;
    }
    
    // Select random instance of this type
    int instanceId = selectRandomInstance(typeId);
    if (instanceId < 0) {
        result.accepted = false;
        // Consistent probability storage: 0 when N=0, -1 when disabled
        result.acceptanceProbability = shouldStoreProbability() ? 0.0 : -1.0;
        return result;
    }
    
    result.residueIndex = instanceId;
    
    // Get instance info before deletion
    FragmentInstance* instance = reservoir_->getInstance(instanceId);
    if (!instance) {
        result.accepted = false;
        return result;
    }
    
    // Save position and orientation for potential restoration
    Vector3 savedPosition = instance->position;
    Quaternion savedOrientation = instance->orientation;
    result.position = savedPosition;

    // Calculate CBMC bias for deletion if enabled
    double cbmcBias = 1.0;
    int numTrials = 1;
    if (useConfBias_ && typeId < static_cast<int>(cbmcTrialsPerType_.size())) {
        numTrials = cbmcTrialsPerType_[typeId];
    }

    if (useConfBias_ && numTrials > 1) {
        // For deletion, calculate W_old: energy of current config + K-1 trial configs
        std::vector<TrialConfiguration> trials;
        trials.reserve(numTrials);

        // Add current configuration as first trial
        TrialConfiguration current;
        current.position = savedPosition;
        current.orientation = savedOrientation;
        current.energy = calculateFragmentEnergy(instanceId);
        trials.push_back(current);

        // Generate K-1 additional trials
        FragmentTemplate* tmpl = reservoir_->getTemplate(typeId);
        if (tmpl) {
            for (int k = 1; k < numTrials; ++k) {
                TrialConfiguration trial;
                trial.position = (useCavityBias_ && cavityManager_) ?
                                generateCavityPosition() : generateRandomPosition();
                applyPeriodicBoundary(trial.position);
                trial.orientation = generateRandomOrientation();

                // Create temporary instance for energy calculation
                int tempId = reservoir_->createInstance(typeId, trial.position, trial.orientation);
                if (tempId >= 0) {
                    synchronizeStateWithReservoir(tempId, true);

                    if (energyMethod_ == EnergyMethod::DIRECT) {
                        cpu::computeResidueEnergyCutoffPBC(*state_, tempId);
                        const auto& residue = state_->residues[tempId];
                        trial.energy = residue.energy_vdw + residue.energy_elec;
                    } else {
                        trial.energy = calculateFragmentEnergy(tempId);
                    }

                    reservoir_->deleteInstance(tempId);
                    synchronizeStateWithReservoir(tempId, false);
                    trials.push_back(trial);
                }
            }
        }

        // Calculate W_old (sum of Boltzmann factors)
        if (trials.size() == static_cast<size_t>(numTrials)) {
            double beta = 1.0 / (8.314e-3 * temperature_);
            double minEnergy = std::numeric_limits<double>::max();
            for (const auto& trial : trials) {
                minEnergy = std::min(minEnergy, trial.energy);
            }

            double sumBoltzmann = 0.0;
            for (const auto& trial : trials) {
                sumBoltzmann += std::exp(-beta * (trial.energy - minEnergy));
            }
            cbmcBias = sumBoltzmann / numTrials * std::exp(beta * minEnergy);  // W_old / K
        }
    }

    // Calculate energy change - optimize for DIRECT mode
    if (energyMethod_ == EnergyMethod::DIRECT) {
        // Fast local ΔE calculation - compute energy of residue to be deleted
        cpu::computeResidueEnergyCutoffPBC(*state_, instanceId);
        const auto& residue = state_->residues[instanceId];
        double residueEnergy = residue.energy_vdw + residue.energy_elec;
        result.deltaE = -residueEnergy;  // Removing this energy from system
        result.energyBefore = residueEnergy;  // For consistency
        result.energyAfter = 0.0;
    } else {
        // Full system energy for Ewald/PME modes
        result.energyBefore = calculateSystemEnergy();
    }
    
    // Temporarily delete (convert to ghost)
    reservoir_->deleteInstance(instanceId);
    synchronizeStateWithReservoir(instanceId, false);
    
    // Calculate energy after deletion for non-DIRECT modes
    if (energyMethod_ != EnergyMethod::DIRECT) {
        result.energyAfter = calculateSystemEnergy();
        result.deltaE = result.energyAfter - result.energyBefore;
    }
    
    // CRITICAL FIX: Calculate deletion bias using saved position for robustness
    // Include CBMC bias in total bias
    result.bias = calculateDeletionBiasAtPosition(savedPosition) * cbmcBias;
    
    // Calculate acceptance probability using proper GCMC formula
    bool accept = false;
    double prob = 0.0;
    if (acceptanceCalculator_) {
        // CRITICAL FIX: Use N_before for correct detailed balance
        prob = acceptanceCalculator_->calculateDeletionProbability(
            typeId, N_before, result.deltaE, result.bias);
        
        // Only store probability if configured
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        
        accept = acceptanceCalculator_->acceptMove(prob);
    } else {
        // Fallback to simple acceptance (should not be used in production)
        double beta = 1.0 / (8.314e-3 * temperature_);
        double activity = 100.0; // Default activity
        prob = std::min(1.0, (N_before / (activity * state_->info.volume)) * 
                       std::exp(-beta * result.deltaE) * result.bias);
        // Apply same probability storage logic as main branch
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        // Use the calculated prob directly for consistency
        accept = (uniform_(rng_) < prob);
    }
    
    if (accept) {
        result.accepted = true;
        acceptedMoves_++;
        energyCache_.invalidate();
    } else {
        // Restore the deleted instance at the same position (keeps same ID)
        bool restored = reservoir_->restoreInstance(instanceId, savedPosition, savedOrientation);
        if (restored) {
            // Successfully restored with same instance ID
            synchronizeStateWithReservoir(instanceId, true);
            result.residueIndex = instanceId;  // Keep original index
        } else {
            // Fallback: create a new instance if restoration failed
            // This can happen if the instance was already purged
            int restoredId = reservoir_->createInstance(typeId, savedPosition, savedOrientation);
            if (restoredId >= 0) {
                synchronizeStateWithReservoir(restoredId, true);
                result.residueIndex = restoredId;
            } else {
                std::cerr << "WARNING: Failed to restore deleted instance in GCMC deletion rejection" << std::endl;
            }
        }
        result.accepted = false;
    }
    
    totalMoves_++;
    
    // Sample statistics if configured using a lightweight countdown
    if (collectStats_) {
        if (--statsCountdown_ <= 0) {
            int particleCount = reservoir_ ? reservoir_->getActiveCount() : 0;
            double energy = calculateSystemEnergy();
            statistics_.addSample(totalMoves_, particleCount, energy,
                                  getAcceptanceRate(), temperature_, 100.0);
            statsCountdown_ = std::max(1, statsInterval_);
        }
    }
    
    return result;
}

// Attempt translation
GCMCEngine::MoveResult GCMCEngine::attemptTranslation(int residueIdx) {
    MoveResult result;
    result.type = MoveResult::TRANSLATE;
    result.residueIndex = residueIdx;
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (!instance || !instance->isActive) {
        result.accepted = false;
        return result;
    }
    
    result.fragmentType = instance->templateId;
    
    // Store old position
    Vector3 oldPos = instance->position;
    result.position = oldPos;
    
    // Calculate energy before move
    result.energyBefore = calculateFragmentEnergy(residueIdx);
    
    // Generate translation using configured step size
    Vector3 displacement = generateTranslationVector(maxTranslationStep_);
    Vector3 newPos = oldPos + displacement;
    applyPeriodicBoundary(newPos);

    // Enforce region constraint: reject moves that leave region
    if (regionConstraint_) {
        movement::Vector3 movNew(newPos.x, newPos.y, newPos.z);
        if (!regionConstraint_->isInRegion(movNew)) {
            // Reject without changing position
            result.deltaE = 0.0;
            result.energyAfter = result.energyBefore;
            result.accepted = false;
            // Probability bookkeeping (only if configured)
            result.acceptanceProbability = shouldStoreProbability() ? 0.0 : -1.0;

            totalMoves_++;
            // Sample statistics if configured using a lightweight countdown
            if (collectStats_) {
                if (--statsCountdown_ <= 0) {
                    int particleCount = reservoir_ ? reservoir_->getActiveCount() : 0;
                    double energy = calculateSystemEnergy();
                    statistics_.addSample(totalMoves_, particleCount, energy,
                                          getAcceptanceRate(), temperature_, 100.0);
                    statsCountdown_ = std::max(1, statsInterval_);
                }
            }
            return result;
        }
    }

    // Update position
    updateFragmentPosition(residueIdx, newPos);
    
    // Calculate energy after move
    result.energyAfter = calculateFragmentEnergy(residueIdx);
    result.deltaE = result.energyAfter - result.energyBefore;
    
    // Accept or reject using unified acceptance calculator
    bool accept = false;
    double prob = 0.0;
    if (acceptanceCalculator_) {
        prob = acceptanceCalculator_->calculateTranslationProbability(result.deltaE, 1.0);
        
        // Only store probability if configured
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        
        accept = acceptanceCalculator_->acceptMove(prob);
    } else {
        // Metropolis criterion
        double beta = 1.0 / (8.314e-3 * temperature_);
        prob = std::min(1.0, std::exp(-beta * result.deltaE));
        // Apply same probability storage logic as main branch
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        accept = acceptMove(result.deltaE, 1.0, temperature_);
    }
    
    if (accept) {
        result.accepted = true;
        result.position = newPos;
        acceptedMoves_++;
        energyCache_.invalidate();
    } else {
        // Restore old position
        updateFragmentPosition(residueIdx, oldPos);
        result.accepted = false;
    }
    
    totalMoves_++;
    
    // Sample statistics if configured using a lightweight countdown
    if (collectStats_) {
        if (--statsCountdown_ <= 0) {
            int particleCount = reservoir_ ? reservoir_->getActiveCount() : 0;
            double energy = calculateSystemEnergy();
            statistics_.addSample(totalMoves_, particleCount, energy,
                                  getAcceptanceRate(), temperature_, 100.0);
            statsCountdown_ = std::max(1, statsInterval_);
        }
    }
    
    return result;
}

// Attempt rotation
GCMCEngine::MoveResult GCMCEngine::attemptRotation(int residueIdx) {
    MoveResult result;
    result.type = MoveResult::ROTATE;
    result.residueIndex = residueIdx;
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (!instance || !instance->isActive) {
        result.accepted = false;
        return result;
    }
    
    result.fragmentType = instance->templateId;
    result.position = instance->position;
    
    // Store old orientation
    Quaternion oldOrient = instance->orientation;
    
    // Calculate energy before rotation
    result.energyBefore = calculateFragmentEnergy(residueIdx);
    
    // Generate rotation using configured angle
    Quaternion rotation = generateRotationQuaternion(maxRotationAngleRad_);
    // CRITICAL FIX: Actually apply the rotation by quaternion multiplication
    Quaternion newOrient = oldOrient * rotation;
    newOrient.normalize();
    
    // Update orientation
    updateFragmentOrientation(residueIdx, newOrient);
    
    // Calculate energy after rotation
    result.energyAfter = calculateFragmentEnergy(residueIdx);
    result.deltaE = result.energyAfter - result.energyBefore;
    
    // Use unified acceptance calculation for rotation
    bool accept = false;
    double prob = 0.0;
    if (acceptanceCalculator_) {
        // Rotation uses standard Metropolis criterion (no N dependence)
        prob = acceptanceCalculator_->calculateTranslationProbability(
            result.deltaE, 1.0);  // bias = 1.0 for rotation
        
        // Only store probability if configured
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        
        accept = acceptanceCalculator_->acceptMove(prob);
    } else {
        // Fallback to direct calculation
        double beta = 1.0 / (8.314e-3 * temperature_);
        prob = (result.deltaE <= 0) ? 1.0 : std::exp(-beta * result.deltaE);
        // Apply same probability storage logic as main branch
        result.acceptanceProbability = shouldStoreProbability() ? prob : -1.0;
        accept = acceptMove(result.deltaE, 1.0, temperature_);
    }
    
    if (accept) {
        result.accepted = true;
        acceptedMoves_++;
        energyCache_.invalidate();
    } else {
        // Restore old orientation
        updateFragmentOrientation(residueIdx, oldOrient);
        result.accepted = false;
    }
    
    totalMoves_++;
    
    // Sample statistics if configured using a lightweight countdown
    if (collectStats_) {
        if (--statsCountdown_ <= 0) {
            int particleCount = reservoir_ ? reservoir_->getActiveCount() : 0;
            double energy = calculateSystemEnergy();
            statistics_.addSample(totalMoves_, particleCount, energy,
                                  getAcceptanceRate(), temperature_, 100.0);
            statsCountdown_ = std::max(1, statsInterval_);
        }
    }
    
    return result;
}

// Attempt swap
GCMCEngine::MoveResult GCMCEngine::attemptSwap(int typeId1, int typeId2) {
    MoveResult result;
    result.type = MoveResult::SWAP;
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    // Select instances to swap
    int idx1 = selectRandomInstance(typeId1);
    int idx2 = selectRandomInstance(typeId2);
    
    if (idx1 < 0 || idx2 < 0) {
        result.accepted = false;
        return result;
    }
    
    // For simplicity, swap is delete type1 + insert type2
    // In practice, would swap identities directly
    
    result.accepted = false;  // Not fully implemented
    totalMoves_++;
    return result;
}

// Attempt regrowth
GCMCEngine::MoveResult GCMCEngine::attemptRegrowth(int residueIdx) {
    MoveResult result;
    result.type = MoveResult::REGROWTH;
    result.residueIndex = residueIdx;
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    // Regrowth = deletion + insertion at new position
    // Simplified implementation
    
    result.accepted = false;  // Not fully implemented
    totalMoves_++;
    return result;
}

// Attempt cluster move
GCMCEngine::MoveResult GCMCEngine::attemptClusterMove(int residueIdx, double cutoff) {
    MoveResult result;
    result.type = MoveResult::CLUSTER;
    result.residueIndex = residueIdx;
    
    if (!state_ || !reservoir_) {
        result.accepted = false;
        return result;
    }
    
    // Find cluster
    std::vector<int> cluster = selectCluster(residueIdx, cutoff);
    
    if (cluster.empty()) {
        result.accepted = false;
        return result;
    }
    
    // Move entire cluster together
    // Simplified - not fully implemented
    
    result.accepted = false;
    totalMoves_++;
    return result;
}

// Select random fragment type
int GCMCEngine::selectRandomFragment() {
    if (!reservoir_) return -1;
    
    int nTypes = reservoir_->getTemplateCount();
    if (nTypes == 0) return -1;
    
    std::uniform_int_distribution<int> dist(0, nTypes - 1);
    return dist(rng_);
}

// Select random instance
int GCMCEngine::selectRandomInstance(int typeId) {
    if (!reservoir_) return -1;

    std::vector<int> instances = reservoir_->getActiveInstances(typeId);
    if (instances.empty()) return -1;

    // If region constraint is set, filter instances to those within the region
    if (regionConstraint_) {
        std::vector<int> filteredInstances;
        for (int id : instances) {
            FragmentInstance* inst = reservoir_->getInstance(id);
            if (inst) {
                movement::Vector3 pos(inst->position.x, inst->position.y, inst->position.z);
                if (regionConstraint_->isInRegion(pos)) {
                    filteredInstances.push_back(id);
                }
            }
        }

        // Use filtered list if not empty, otherwise fall back to all instances
        // (this prevents deletion from being completely blocked if all molecules drift outside)
        if (!filteredInstances.empty()) {
            instances = filteredInstances;
        }
    }

    std::uniform_int_distribution<int> dist(0, instances.size() - 1);
    return instances[dist(rng_)];
}

// Select cluster
std::vector<int> GCMCEngine::selectCluster(int seedIdx, double cutoff) {
    std::vector<int> cluster;
    if (!reservoir_) return cluster;
    
    FragmentInstance* seed = reservoir_->getInstance(seedIdx);
    if (!seed) return cluster;
    
    cluster.push_back(seedIdx);
    
    // Find neighbors within cutoff
    std::vector<int> allInstances = reservoir_->getActiveInstances();
    for (int idx : allInstances) {
        if (idx == seedIdx) continue;
        
        FragmentInstance* instance = reservoir_->getInstance(idx);
        if (!instance) continue;
        
        double dist = minimumImageDistance(seed->position, instance->position);
        if (dist < cutoff) {
            cluster.push_back(idx);
        }
    }
    
    return cluster;
}

// Generate random position
Vector3 GCMCEngine::generateRandomPosition() {
    if (!state_) {
        return Vector3(uniform_(rng_) * 100,
                      uniform_(rng_) * 100,
                      uniform_(rng_) * 100);
    }
    
    // Get box dimensions with fallback
    double boxX, boxY, boxZ;
    if (state_->periodicBox.size() >= 3) {
        boxX = state_->periodicBox[0];
        boxY = state_->periodicBox[1];
        boxZ = state_->periodicBox[2];
    } else if (state_->info.box[0] > 0 && state_->info.box[1] > 0 && state_->info.box[2] > 0) {
        // Fallback to info.box if periodicBox not set
        boxX = state_->info.box[0];
        boxY = state_->info.box[1];
        boxZ = state_->info.box[2];
    } else {
        // No valid box dimensions, use default
        boxX = boxY = boxZ = 100.0;
    }
    
    // If region constraint is set, sample from constrained region
    if (regionConstraint_) {
        movement::Vector3 movPos = regionConstraint_->samplePosition();
        return Vector3(movPos.x, movPos.y, movPos.z);
    }

    // Use [0, L) coordinate system to match CavityManager
    return Vector3(uniform_(rng_) * boxX,
                  uniform_(rng_) * boxY,
                  uniform_(rng_) * boxZ);
}

// Generate cavity position
Vector3 GCMCEngine::generateCavityPosition() {
    if (!cavityManager_) {
        return generateRandomPosition();
    }

    // Try to find a cavity position within the region constraint
    const int maxAttempts = 100;
    for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        // Get cavity from manager (returns movement::Vector3)
        movement::Vector3 movCavity = cavityManager_->selectCavity();

        // Convert to montecarlo::Vector3
        Vector3 cavity(movCavity.x, movCavity.y, movCavity.z);

        // Add small random displacement
        Vector3 displacement(normal_(rng_) * 0.5,
                            normal_(rng_) * 0.5,
                            normal_(rng_) * 0.5);

        Vector3 position = cavity + displacement;

        // Check if position is within region constraint
        if (!regionConstraint_ || regionConstraint_->isInRegion(movement::Vector3(position.x, position.y, position.z))) {
            return position;
        }
    }

    // If no valid cavity found in region, fall back to random position in region
    if (regionConstraint_) {
        movement::Vector3 movPos = regionConstraint_->samplePosition();
        return Vector3(movPos.x, movPos.y, movPos.z);
    }

    // Last resort: random position in box
    return generateRandomPosition();
}

// Generate random orientation
Quaternion GCMCEngine::generateRandomOrientation() {
    // Use the correct Shoemake algorithm for uniform quaternion distribution
    // This matches gcmc_gpu's create_random_quarternion() implementation
    double u = uniform_(rng_);
    double v = uniform_(rng_);
    double w = uniform_(rng_);
    
    double sqrt_1_minus_u = std::sqrt(1.0 - u);
    double sqrt_u = std::sqrt(u);
    double two_pi_v = 2.0 * M_PI * v;
    double two_pi_w = 2.0 * M_PI * w;
    
    Quaternion q(
        sqrt_u * std::cos(two_pi_w),           // w component
        sqrt_1_minus_u * std::sin(two_pi_v),   // x component
        sqrt_1_minus_u * std::cos(two_pi_v),   // y component
        sqrt_u * std::sin(two_pi_w)            // z component
    );
    q.normalize();
    
    return q;
}

// Generate translation vector
Vector3 GCMCEngine::generateTranslationVector(double maxDist) {
    // Random direction
    double theta = uniform_(rng_) * 2 * M_PI;
    double phi = std::acos(2 * uniform_(rng_) - 1);
    
    // Random magnitude
    double r = uniform_(rng_) * maxDist;
    
    return Vector3(r * std::sin(phi) * std::cos(theta),
                  r * std::sin(phi) * std::sin(theta),
                  r * std::cos(phi));
}

// Generate rotation quaternion
Quaternion GCMCEngine::generateRotationQuaternion(double maxAngle) {
    // Random axis
    Vector3 axis(normal_(rng_), normal_(rng_), normal_(rng_));
    double norm = axis.norm();
    if (norm > 0) {
        axis = axis * (1.0 / norm);
    }
    
    // Random angle
    double angle = uniform_(rng_) * maxAngle;
    
    // Create quaternion from axis-angle
    double halfAngle = angle / 2;
    double s = std::sin(halfAngle);
    
    Quaternion q(std::cos(halfAngle), s * axis.x, s * axis.y, s * axis.z);
    q.normalize();
    
    return q;
}

// Calculate system energy
double GCMCEngine::calculateSystemEnergy() {
    if (!state_) return 0.0;
    
    if (energyCache_.valid) {
        return energyCache_.totalEnergy;
    }
    
    double totalEnergy = 0.0;
    
    // Use energy callback if available
    if (energyCallback_) {
        totalEnergy = energyCallback_->calculateSystemEnergy(*state_);
    } else {
        // Fallback to direct energy module usage
        if (energyMethod_ == EnergyMethod::PME) {
            computeSystemEnergy(*state_, EnergyMethod::PME);
        } else if (energyMethod_ == EnergyMethod::EWALD) {
            computeSystemEnergy(*state_, EnergyMethod::EWALD);
        } else {
            // DIRECT with cutoff and PBC
            computeSystemEnergy(*state_, EnergyMethod::DIRECT, true, true);
        }
        
        // Get total energy from state
        if (energyMethod_ == EnergyMethod::EWALD || energyMethod_ == EnergyMethod::PME) {
            totalEnergy = energy::getTotalEnergy(*state_, energyMethod_);
        } else {
            totalEnergy = energy::getTotalEnergy(*state_, EnergyMethod::DIRECT);
        }
    }
    
    energyCache_.totalEnergy = totalEnergy;
    energyCache_.valid = true;
    
    return totalEnergy;
}

// Calculate fragment energy
double GCMCEngine::calculateFragmentEnergy(int residueIdx) {
    if (!state_ || residueIdx < 0 || residueIdx >= static_cast<int>(state_->residues.size())) {
        return 0.0;
    }
    
    // Get residue from state
    auto& residue = state_->residues[residueIdx];
    if (!residue.active) return 0.0;
    
    // Check cache
    if (reservoir_) {
        FragmentInstance* instance = reservoir_->getInstance(residueIdx);
        if (instance && instance->lastEnergyUpdate == totalMoves_) {
            return instance->energy_total;
        }
    }
    
    double energy = 0.0;
    
    // Use energy callback if available
    if (energyCallback_) {
        energy = energyCallback_->calculateResidueEnergy(*state_, residueIdx);
    } else {
        // Fallback to direct energy module usage
        // Mark residue as moved for movement energy calculation
        // residue.moved /* moved flag not in MCResidue */ = true;
        
        // Calculate movement energy for this residue
        if (energyMethod_ == EnergyMethod::PME) {
            computeMovementEnergy(*state_, EnergyMethod::PME);
        } else if (energyMethod_ == EnergyMethod::EWALD) {
            computeMovementEnergy(*state_, EnergyMethod::EWALD);
        } else {
            computeMovementEnergy(*state_, EnergyMethod::DIRECT, true, true);
        }
        
        // Reset moved flag
        // residue.moved /* moved flag not in MCResidue */ = false;
        
        energy = residue.energy_vdw + residue.energy_elec;
    }
    
    // Update cache
    if (reservoir_) {
        FragmentInstance* instance = reservoir_->getInstance(residueIdx);
        if (instance) {
            instance->energy_total = energy;
            instance->lastEnergyUpdate = totalMoves_;
        }
    }
    
    return energy;
}

// Calculate interaction energy
double GCMCEngine::calculateInteractionEnergy(int residueIdx) {
    // This function is now deprecated - use calculateFragmentEnergy instead
    // which properly uses the energy module
    return calculateFragmentEnergy(residueIdx);
}

// Calculate pair energy
double GCMCEngine::calculatePairEnergy(int idx1, int idx2) {
    // Suppress unused parameter warnings
    (void)idx1;
    (void)idx2;
    
    // This function is now deprecated - the energy module handles
    // all pair interactions properly with the correct force field parameters
    // Use calculateFragmentEnergy or calculateSystemEnergy instead
    return 0.0;
}

// Calculate insertion bias
double GCMCEngine::calculateInsertionBias(const FragmentTemplate& tmpl,
                                         const Vector3& position,
                                         const Quaternion& orientation) {
    // Suppress unused parameter warnings
    (void)tmpl;
    (void)orientation;

    double bias = 1.0;

    // CRITICAL: Only apply cavity bias if enabled
    if (cavityManager_ && useCavityBias_) {
        // Use position-based cavity score (O(1)) instead of volume calculation (O(n³))
        // This is much faster and was the original implementation
        movement::Vector3 pos(position.x, position.y, position.z);
        bias *= cavityManager_->getCavityScore(pos);
    }

    if (configBias_) {
        // Config bias calculation would go here
        bias *= 1.0;
    }

    // Proposal bias for detailed balance when using target_numwaters
    // For insertion: multiply by p_delete/p_insert ratio
    double proposalBias = getConfigValue("proposalBias");
    if (proposalBias > 0) {
        bias *= proposalBias;  // This is p_delete/p_insert for insertion
    }

    return bias;
}

// Calculate deletion bias at specific position (more robust)
double GCMCEngine::calculateDeletionBiasAtPosition(const Vector3& position) {
    double bias = 1.0;

    // CRITICAL: For detailed balance, deletion bias must match insertion bias
    // at the same position
    if (cavityManager_ && useCavityBias_) {
        movement::Vector3 pos(position.x, position.y, position.z);
        // Use the same cavity score calculation as insertion
        bias *= cavityManager_->getCavityScore(pos);
    }

    if (configBias_) {
        // Config bias calculation would go here (must match insertion)
        bias *= 1.0;
    }

    // Proposal bias for detailed balance when using target_numwaters
    // For deletion: multiply by p_insert/p_delete ratio (inverse of insertion)
    double proposalBias = getConfigValue("proposalBias");
    if (proposalBias > 0) {
        bias *= (1.0 / proposalBias);  // This is p_insert/p_delete for deletion
    }

    return bias;
}

// Calculate deletion bias (legacy - depends on reservoir state)
double GCMCEngine::calculateDeletionBias(int residueIdx) {
    // Try to get position from reservoir
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (instance) {
        Vector3 pos(instance->position.x, instance->position.y, instance->position.z);
        return calculateDeletionBiasAtPosition(pos);
    }
    // Fallback if instance not accessible
    return 1.0;
}

// Calculate regrowth bias
double GCMCEngine::calculateRegrowthBias(int residueIdx) {
    // Combination of deletion and insertion biases
    double deletionBias = calculateDeletionBias(residueIdx);
    
    // Would calculate insertion bias at new position
    double insertionBias = 1.0;
    
    return deletionBias * insertionBias;
}

// Accept move
bool GCMCEngine::acceptMove(double deltaE, double bias, double temperature) {
    if (deltaE <= 0) return true;
    
    double beta = 1.0 / (8.314e-3 * temperature);  // kJ/(mol*K)
    double probability = bias * std::exp(-beta * deltaE);
    
    return uniform_(rng_) < probability;
}

// Calculate acceptance probability
double GCMCEngine::calculateAcceptanceProbability(const MoveResult& result,
                                                 double temperature) {
    double beta = 1.0 / (8.314e-3 * temperature);
    return std::min(1.0, result.bias * std::exp(-beta * result.deltaE));
}

// Synchronize MCState with reservoir
void GCMCEngine::synchronizeStateWithReservoir(int instanceId, bool isInsertion) {
    if (!state_ || !reservoir_) return;
    
    FragmentInstance* instance = reservoir_->getInstance(instanceId);
    if (!instance) return;
    
    if (isInsertion) {
        // Add residue to MCState
        if (instanceId >= static_cast<int>(state_->residues.size())) {
            // Need to expand residues vector
            state_->residues.resize(instanceId + 1);
        }
        
        // Get template for atom information
        const FragmentTemplate* tmpl = reservoir_->getTemplate(instance->templateId);
        if (!tmpl) return;
        
        // Update residue in state
        auto& residue = state_->residues[instanceId];
        residue.active = true;
        residue.resid = instanceId;
        residue.resname = tmpl->name;
        // MCResidue doesn't have chainid field
        // residue.moved /* moved flag not in MCResidue */ = false;
        residue.energy_vdw = 0.0;
        residue.energy_elec = 0.0;
        
        // CRITICAL: Set atomStart and atomCount for energy calculations
        residue.atomStart = state_->activeAtomCount;
        residue.atomCount = static_cast<int>(tmpl->atoms.size());
        
        // Clear and add atoms
        residue.atoms.clear();
        residue.atoms.reserve(tmpl->atoms.size());
        
        // Transform template atoms by instance position and orientation
        for (const auto& tmplAtom : tmpl->atoms) {
            MCAtom atom;
            // MCAtom doesn't have active or atomId fields
            // atom.active = true;
            // atom.atomId = tmplAtom.id;  // tmplAtom may not have id field
            atom.type = tmplAtom.type;
            atom.charge = tmplAtom.charge;
            atom.mass = tmplAtom.mass;
            atom.name = tmplAtom.name;
            
            // Apply rotation and translation
            // Quaternion rotation of vector: q * v * q^-1
            // For unit quaternion, q^-1 = (w, -x, -y, -z)
            Vector3 v(tmplAtom.x, tmplAtom.y, tmplAtom.z);
            Quaternion q = instance->orientation;
            
            // Use the movement layer's Quaternion::rotate() method
            double vx = v.x, vy = v.y, vz = v.z;
            Vector3 templatePos(vx, vy, vz);
            Vector3 rotatedPos = q.rotate(templatePos);
            
            atom.x = instance->position.x + rotatedPos.x;
            atom.y = instance->position.y + rotatedPos.y;
            atom.z = instance->position.z + rotatedPos.z;
            
            // Sync position vector with x,y,z coordinates
            atom.updatePosition();
            
            residue.atoms.push_back(atom);
            
            // CRITICAL: Also add to global atoms array for energy calculations
            if (state_->activeAtomCount < static_cast<int>(state_->atoms.size())) {
                state_->atoms[state_->activeAtomCount] = atom;
            } else {
                state_->atoms.push_back(atom);
            }
            state_->activeAtomCount++;
        }
        
        // Update residue index in fragment instance
        instance->residueIndex = instanceId;
        
        // CRITICAL: Update activeResidueCount
        // Find the highest active residue index + 1
        state_->activeResidueCount = 0;
        for (int i = 0; i < static_cast<int>(state_->residues.size()); i++) {
            if (state_->residues[i].active) {
                state_->activeResidueCount = i + 1;
            }
        }
    } else {
        // Mark residue as inactive (deletion case)
        if (instanceId < static_cast<int>(state_->residues.size())) {
            auto& residue = state_->residues[instanceId];
            residue.active = false;
            
            // CRITICAL: Don't actually remove atoms from global array to avoid shifting indices
            // Just mark the residue as inactive so energy calculations skip it
            residue.atomCount = 0;  // Mark as having no atoms
            residue.atoms.clear();
            
            // CRITICAL: Update activeResidueCount
            // Find the highest active residue index + 1
            state_->activeResidueCount = 0;
            for (int i = 0; i < static_cast<int>(state_->residues.size()); i++) {
                if (state_->residues[i].active) {
                    state_->activeResidueCount = i + 1;
                }
            }
            
            // Note: We don't decrement activeAtomCount here to avoid index shifting
            // This is a simplification for now - a production system would compact arrays
        }
    }
}

// Update fragment position
void GCMCEngine::updateFragmentPosition(int residueIdx, const Vector3& newPos) {
    if (!reservoir_) return;
    
    reservoir_->updatePosition(residueIdx, newPos);
    energyCache_.invalidate();
    
    // Update atom coordinates without adding new atoms
    updateAtomCoordinates(residueIdx);
}

// Update fragment orientation
void GCMCEngine::updateFragmentOrientation(int residueIdx, const Quaternion& newOrient) {
    if (!reservoir_) return;
    
    reservoir_->updateOrientation(residueIdx, newOrient);
    energyCache_.invalidate();
    
    // Update atom coordinates without adding new atoms
    updateAtomCoordinates(residueIdx);
}

// Update atom coordinates for an existing residue without changing atom count
void GCMCEngine::updateAtomCoordinates(int residueIdx) {
    if (!state_ || !reservoir_) return;
    
    if (residueIdx >= static_cast<int>(state_->residues.size())) return;
    
    auto& residue = state_->residues[residueIdx];
    if (!residue.active) return;
    
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (!instance) return;
    
    const FragmentTemplate* tmpl = reservoir_->getTemplate(instance->templateId);
    if (!tmpl) return;
    
    // Update atoms in both residue.atoms and state->atoms arrays
    int atomIdx = 0;
    for (const auto& tmplAtom : tmpl->atoms) {
        if (atomIdx >= residue.atomCount) break;
        
        // Apply rotation and translation
        Vector3 v(tmplAtom.x, tmplAtom.y, tmplAtom.z);
        Quaternion q = instance->orientation;
        
        // Use the verified Quaternion::rotate() method
        Vector3 rotatedPos = q.rotate(v);
        
        // Update coordinates in residue.atoms
        if (atomIdx < static_cast<int>(residue.atoms.size())) {
            residue.atoms[atomIdx].x = instance->position.x + rotatedPos.x;
            residue.atoms[atomIdx].y = instance->position.y + rotatedPos.y;
            residue.atoms[atomIdx].z = instance->position.z + rotatedPos.z;
            residue.atoms[atomIdx].updatePosition();
        }
        
        // Update coordinates in global atoms array
        int globalIdx = residue.atomStart + atomIdx;
        if (globalIdx < static_cast<int>(state_->atoms.size())) {
            state_->atoms[globalIdx].x = instance->position.x + rotatedPos.x;
            state_->atoms[globalIdx].y = instance->position.y + rotatedPos.y;
            state_->atoms[globalIdx].z = instance->position.z + rotatedPos.z;
            state_->atoms[globalIdx].updatePosition();
        }
        
        atomIdx++;
    }
}

// Apply periodic boundary conditions
void GCMCEngine::applyPeriodicBoundary(Vector3& position) {
    if (!state_) return;
    
    // Get box dimensions with fallback
    double boxX, boxY, boxZ;
    if (state_->periodicBox.size() >= 3) {
        boxX = state_->periodicBox[0];
        boxY = state_->periodicBox[1];
        boxZ = state_->periodicBox[2];
    } else if (state_->info.box[0] > 0 && state_->info.box[1] > 0 && state_->info.box[2] > 0) {
        // Fallback to info.box if periodicBox not set
        boxX = state_->info.box[0];
        boxY = state_->info.box[1];
        boxZ = state_->info.box[2];
    } else {
        // No valid box dimensions, skip PBC
        return;
    }
    
    // Use [0, L) coordinate system
    while (position.x < 0) position.x += boxX;
    while (position.x >= boxX) position.x -= boxX;
    while (position.y < 0) position.y += boxY;
    while (position.y >= boxY) position.y -= boxY;
    while (position.z < 0) position.z += boxZ;
    while (position.z >= boxZ) position.z -= boxZ;
}

// Calculate minimum image distance
double GCMCEngine::minimumImageDistance(const Vector3& r1, const Vector3& r2) {
    if (!state_) {
        return (r1 - r2).norm();
    }
    
    Vector3 dr = r1 - r2;
    
    // Get box dimensions with fallback (same as applyPeriodicBoundary)
    double boxX, boxY, boxZ;
    if (state_->periodicBox.size() >= 3) {
        boxX = state_->periodicBox[0];
        boxY = state_->periodicBox[1];
        boxZ = state_->periodicBox[2];
    } else if (state_->info.box[0] > 0 && state_->info.box[1] > 0 && state_->info.box[2] > 0) {
        // Fallback to info.box if periodicBox not set
        boxX = state_->info.box[0];
        boxY = state_->info.box[1];
        boxZ = state_->info.box[2];
    } else {
        // No valid box dimensions, return direct distance
        return dr.norm();
    }
    
    // Apply minimum image convention for [0, L) coordinate system
    if (std::abs(dr.x) > boxX/2) {
        dr.x = dr.x - std::copysign(boxX, dr.x);
    }
    if (std::abs(dr.y) > boxY/2) {
        dr.y = dr.y - std::copysign(boxY, dr.y);
    }
    if (std::abs(dr.z) > boxZ/2) {
        dr.z = dr.z - std::copysign(boxZ, dr.z);
    }
    
    return dr.norm();
}

// Dynamic configuration implementation
void GCMCEngine::setConfigValue(const std::string& key, double value) {
    configMap_[key] = value;
    
    // Apply specific configuration changes
    if (key == "temperature") {
        temperature_ = value;
        if (acceptanceCalculator_) {
            acceptanceCalculator_->setTemperature(value);
        }
    } else if (key == "cutoff") {
        cutoff_ = value;
    } else if (key == "statsInterval") {
        statsInterval_ = static_cast<int>(value);
        statistics_.setSamplingInterval(statsInterval_);
    } else if (key == "autoAdjustStats") {
        statistics_.setAutoAdjust(value > 0.5);
    } else if (key == "collectStats") {
        collectStats_ = (value > 0.5);
    } else if (key == "maxTranslation") {
        maxTranslationStep_ = value;
    } else if (key == "maxRotation" || key == "maxRotationAngle") {
        maxRotationAngleRad_ = value;  // Support both key names
    } else if (key == "useCavityBias") {
        useCavityBias_ = (value > 0.5);
    } else if (key == "useConfBias") {
        useConfBias_ = (value > 0.5);
    } else if (key == "storeProbabilities") {
        // Clear cache when config changes
        storeProbabilityCached_ = false;
    } else if (key == "storeAcceptanceProbability") {
        // Clear cache when config changes
        storeProbabilityCached_ = false;
    }
}

double GCMCEngine::getConfigValue(const std::string& key) const {
    auto it = configMap_.find(key);
    if (it != configMap_.end()) {
        return it->second;
    }
    
    // Return current values for known keys
    if (key == "temperature") return temperature_;
    if (key == "cutoff") return cutoff_;
    if (key == "statsInterval") return static_cast<double>(statsInterval_);
    if (key == "collectStats") return collectStats_ ? 1.0 : 0.0;
    if (key == "maxTranslation") return maxTranslationStep_;
    if (key == "maxRotation" || key == "maxRotationAngle") return maxRotationAngleRad_;
    if (key == "useCavityBias") return useCavityBias_ ? 1.0 : 0.0;
    if (key == "useConfBias") return useConfBias_ ? 1.0 : 0.0;
    if (key == "storeProbabilities") return shouldStoreProbability() ? 1.0 : 0.0;
    
    return 0.0;  // Default for unknown keys
}

void GCMCEngine::setStatisticsInterval(int interval) {
    statsInterval_ = std::max(1, interval);
    statistics_.setSamplingInterval(statsInterval_);
    configMap_["statsInterval"] = static_cast<double>(statsInterval_);
    // Reset countdown to align with the new interval
    statsCountdown_ = statsInterval_;
}

// Check if should store probability - controlled by configuration
bool GCMCEngine::shouldStoreProbability() const {
    // Use cached value for performance
    if (!storeProbabilityCached_) {
        // First time: check environment and config directly
        // Avoid calling the global function to prevent cross-test contamination
        const char* env_value = std::getenv("GCMC_STORE_PROB");
        if (env_value != nullptr) {
            storeProbabilityValue_ = true;
        } else {
            // Check runtime config keys
            auto it = configMap_.find("storeProbabilities");
            if (it != configMap_.end()) {
                storeProbabilityValue_ = (it->second > 0.5);
            } else {
                auto it2 = configMap_.find("storeAcceptanceProbability");
                if (it2 != configMap_.end()) {
                    storeProbabilityValue_ = (it2->second > 0.5);
                } else {
                    // Default: don't store probabilities for performance
                    storeProbabilityValue_ = false;
                }
            }
        }
        storeProbabilityCached_ = true;
    }
    return storeProbabilityValue_;
}

// CBMC insertion - generate K trials and select based on Boltzmann weights
GCMCEngine::TrialConfiguration GCMCEngine::performCBMCInsertion(int typeId, int numTrials) {
    std::vector<TrialConfiguration> trials;
    trials.reserve(numTrials);

    FragmentTemplate* tmpl = reservoir_->getTemplate(typeId);
    if (!tmpl) {
        // Return default configuration if template not found
        return TrialConfiguration{Vector3(0,0,0), Quaternion(1,0,0,0), 0.0, 1.0};
    }

    // Generate K trial configurations
    double minEnergy = std::numeric_limits<double>::max();
    for (int k = 0; k < numTrials; ++k) {
        TrialConfiguration trial;

        // Generate position and orientation
        trial.position = (useCavityBias_ && cavityManager_) ?
                        generateCavityPosition() : generateRandomPosition();
        applyPeriodicBoundary(trial.position);
        trial.orientation = generateRandomOrientation();

        // Create temporary instance for energy calculation
        int tempId = reservoir_->createInstance(typeId, trial.position, trial.orientation);
        if (tempId < 0) continue;

        synchronizeStateWithReservoir(tempId, true);

        // Calculate energy for this configuration
        if (energyMethod_ == EnergyMethod::DIRECT) {
            cpu::computeResidueEnergyCutoffPBC(*state_, tempId);
            const auto& residue = state_->residues[tempId];
            trial.energy = residue.energy_vdw + residue.energy_elec;
        } else {
            trial.energy = calculateFragmentEnergy(tempId);
        }

        // Clean up temporary instance
        reservoir_->deleteInstance(tempId);
        synchronizeStateWithReservoir(tempId, false);

        // Track minimum energy for numerical stability
        minEnergy = std::min(minEnergy, trial.energy);
        trials.push_back(trial);
    }

    // Calculate Boltzmann weights (subtract minEnergy for numerical stability)
    double beta = 1.0 / (8.314e-3 * temperature_);
    double totalWeight = 0.0;
    for (auto& trial : trials) {
        trial.weight = std::exp(-beta * (trial.energy - minEnergy));
        totalWeight += trial.weight;
    }

    // Select configuration based on weights
    double r = uniform_(rng_) * totalWeight;
    double cumWeight = 0.0;
    for (const auto& trial : trials) {
        cumWeight += trial.weight;
        if (cumWeight >= r) {
            return trial;
        }
    }

    // Fallback to last trial (should not happen)
    return trials.back();
}

// Calculate CBMC bias factor
double GCMCEngine::calculateCBMCBias(const std::vector<TrialConfiguration>& trials, int selectedIdx) {
    if (trials.empty() || selectedIdx < 0 || selectedIdx >= static_cast<int>(trials.size())) {
        return 1.0;
    }

    // Calculate average Boltzmann factor
    double beta = 1.0 / (8.314e-3 * temperature_);
    double minEnergy = std::numeric_limits<double>::max();
    for (const auto& trial : trials) {
        minEnergy = std::min(minEnergy, trial.energy);
    }

    double sumBoltzmann = 0.0;
    for (const auto& trial : trials) {
        sumBoltzmann += std::exp(-beta * (trial.energy - minEnergy));
    }

    // Return W_new / K for insertion
    return sumBoltzmann / trials.size() * std::exp(beta * minEnergy);
}

// Get residue position
Vector3 GCMCEngine::getResiduePosition(int residueIdx) {
    if (!reservoir_) return Vector3(0, 0, 0);
    
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (!instance) return Vector3(0, 0, 0);
    
    return Vector3(instance->position.x, instance->position.y, instance->position.z);
}

// Get residue orientation
Quaternion GCMCEngine::getResidueOrientation(int residueIdx) {
    if (!reservoir_) return Quaternion(1, 0, 0, 0);
    
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (!instance) return Quaternion(1, 0, 0, 0);
    
    return instance->orientation;
}

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc
