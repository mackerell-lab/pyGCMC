#include "GCMCEngine.hpp"
#include "GCMCAcceptance.hpp"
#include "GCMCConfig.hpp"
#include "../../energy/EnergyModule.hpp"
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
    
    // Generate position and orientation based on configuration
    Vector3 position = (useCavityBias_ && cavityManager_) ? 
                      generateCavityPosition() : generateRandomPosition();
    applyPeriodicBoundary(position);  // Ensure position is within PBC
    Quaternion orientation = generateRandomOrientation();
    
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
    
    // Calculate energy after insertion
    result.energyAfter = calculateSystemEnergy();
    result.deltaE = result.energyAfter - result.energyBefore;
    
    // Calculate bias
    result.bias = calculateInsertionBias(*tmpl, position, orientation);
    
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
    
    // Calculate energy before deletion
    result.energyBefore = calculateSystemEnergy();
    
    // Temporarily delete (convert to ghost)
    reservoir_->deleteInstance(instanceId);
    synchronizeStateWithReservoir(instanceId, false);
    
    // Calculate energy after deletion
    result.energyAfter = calculateSystemEnergy();
    result.deltaE = result.energyAfter - result.energyBefore;
    
    // CRITICAL FIX: Calculate deletion bias AFTER deletion for cavity bias
    result.bias = calculateDeletionBias(instanceId);
    
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
    
    // Get cavity from manager (returns movement::Vector3)
    movement::Vector3 movCavity = cavityManager_->selectCavity();
    
    // Convert to montecarlo::Vector3
    Vector3 cavity(movCavity.x, movCavity.y, movCavity.z);
    
    // Add small random displacement
    Vector3 displacement(normal_(rng_) * 0.5,
                        normal_(rng_) * 0.5,
                        normal_(rng_) * 0.5);
    
    return cavity + displacement;
}

// Generate random orientation
Quaternion GCMCEngine::generateRandomOrientation() {
    double u1 = uniform_(rng_);
    double u2 = uniform_(rng_);
    double u3 = uniform_(rng_);
    
    Quaternion q(
        std::sqrt(1 - u1) * std::sin(2 * M_PI * u2),
        std::sqrt(1 - u1) * std::cos(2 * M_PI * u2),
        std::sqrt(u1) * std::sin(2 * M_PI * u3),
        std::sqrt(u1) * std::cos(2 * M_PI * u3)
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
    (void)position;  // Not used for volume-based bias
    
    double bias = 1.0;
    
    if (cavityManager_) {
        // Use position-based cavity score (O(1)) instead of volume calculation (O(n³))
        // This is much faster and was the original implementation
        movement::Vector3 pos(position.x, position.y, position.z);
        bias *= cavityManager_->getCavityScore(pos);
    }
    
    if (configBias_) {
        // Config bias calculation would go here
        bias *= 1.0;
    }
    
    return bias;
}

// Calculate deletion bias
double GCMCEngine::calculateDeletionBias(int residueIdx) {
    // Suppress unused parameter warning
    (void)residueIdx;
    
    double bias = 1.0;
    
    // For deletion, the bias is typically the inverse of insertion bias
    // But since we use position-based scoring, we keep it simple
    // The detailed balance is maintained by the acceptance probability calculation
    
    return bias;
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
            
            // Simplified rotation formula for unit quaternion
            double qw = q.w, qx = q.x, qy = q.y, qz = q.z;
            double vx = v.x, vy = v.y, vz = v.z;
            
            // Calculate rotated position using quaternion rotation formula
            double rx = vx * (qw*qw + qx*qx - qy*qy - qz*qz) + 
                       vy * 2*(qx*qy - qw*qz) + 
                       vz * 2*(qx*qz + qw*qy);
            double ry = vx * 2*(qx*qy + qw*qz) + 
                       vy * (qw*qw - qx*qx + qy*qy - qz*qz) + 
                       vz * 2*(qy*qz - qw*qx);
            double rz = vx * 2*(qx*qz - qw*qy) + 
                       vy * 2*(qy*qz + qw*qx) + 
                       vz * (qw*qw - qx*qx - qy*qy + qz*qz);
            
            atom.x = instance->position.x + rx;
            atom.y = instance->position.y + ry;
            atom.z = instance->position.z + rz;
            
            // Sync position vector with x,y,z coordinates
            atom.updatePosition();
            
            residue.atoms.push_back(atom);
        }
        
        // Update residue index in fragment instance
        instance->residueIndex = instanceId;
    } else {
        // Mark residue as inactive
        if (instanceId < static_cast<int>(state_->residues.size())) {
            state_->residues[instanceId].active = false;
            state_->residues[instanceId].atoms.clear();
        }
    }
}

// Update fragment position
void GCMCEngine::updateFragmentPosition(int residueIdx, const Vector3& newPos) {
    if (!reservoir_) return;
    
    reservoir_->updatePosition(residueIdx, newPos);
    energyCache_.invalidate();
    
    // CRITICAL: Force synchronize MCState after position update
    // This ensures energy calculations and cavity bias see updated positions
    synchronizeStateWithReservoir(residueIdx, true);
}

// Update fragment orientation
void GCMCEngine::updateFragmentOrientation(int residueIdx, const Quaternion& newOrient) {
    if (!reservoir_) return;
    
    reservoir_->updateOrientation(residueIdx, newOrient);
    energyCache_.invalidate();
    
    // CRITICAL: Force synchronize MCState after orientation update
    // This ensures energy calculations see correctly oriented atoms
    synchronizeStateWithReservoir(residueIdx, true);
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
    } else if (key == "maxRotation") {
        maxRotationAngleRad_ = value;
    } else if (key == "useCavityBias") {
        useCavityBias_ = (value > 0.5);
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
    if (key == "maxRotation") return maxRotationAngleRad_;
    if (key == "useCavityBias") return useCavityBias_ ? 1.0 : 0.0;
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
