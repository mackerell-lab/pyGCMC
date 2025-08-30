#include "GCMCEngine.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {
namespace gcmc {

// Constructor
GCMCEngine::GCMCEngine()
    : state_(nullptr),
      reservoir_(nullptr),
      cavityManager_(nullptr),
      configBias_(nullptr),
      temperature_(300.0),
      cutoff_(12.0),
      energyMethod_(EnergyMethod::DIRECT),
      rng_(std::random_device{}()),
      uniform_(0.0, 1.0),
      normal_(0.0, 1.0),
      totalMoves_(0),
      acceptedMoves_(0) {
    energyCache_.valid = false;
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

// Set seed
void GCMCEngine::setSeed(unsigned int seed) {
    rng_.seed(seed);
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
    
    // Generate position and orientation
    Vector3 position = cavityManager_ ? generateCavityPosition() : generateRandomPosition();
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
    
    // Calculate energy after insertion
    result.energyAfter = calculateSystemEnergy();
    result.deltaE = result.energyAfter - result.energyBefore;
    
    // Calculate bias
    result.bias = calculateInsertionBias(*tmpl, position, orientation);
    
    // Accept or reject
    bool accept = acceptMove(result.deltaE, result.bias, temperature_);
    
    if (accept) {
        result.accepted = true;
        result.residueIndex = instanceId;
        acceptedMoves_++;
        energyCache_.invalidate();
    } else {
        // Remove the instance
        reservoir_->deleteInstance(instanceId);
        result.accepted = false;
    }
    
    totalMoves_++;
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
    
    // Select random instance of this type
    int instanceId = selectRandomInstance(typeId);
    if (instanceId < 0) {
        result.accepted = false;
        return result;
    }
    
    result.residueIndex = instanceId;
    
    // Get instance info before deletion
    FragmentInstance* instance = reservoir_->getInstance(instanceId);
    if (!instance) {
        result.accepted = false;
        return result;
    }
    
    result.position = instance->position;
    
    // Calculate energy before deletion
    result.energyBefore = calculateSystemEnergy();
    
    // Calculate deletion bias
    result.bias = calculateDeletionBias(instanceId);
    
    // Temporarily delete (convert to ghost)
    reservoir_->deleteInstance(instanceId);
    
    // Calculate energy after deletion
    result.energyAfter = calculateSystemEnergy();
    result.deltaE = result.energyAfter - result.energyBefore;
    
    // Accept or reject
    bool accept = acceptMove(result.deltaE, result.bias, temperature_);
    
    if (accept) {
        result.accepted = true;
        acceptedMoves_++;
        energyCache_.invalidate();
    } else {
        // Restore the instance (this would need implementation in FragmentReservoir)
        // For now, we'll create a new instance at the same position
        reservoir_->createInstance(typeId, result.position, Quaternion());
        result.accepted = false;
    }
    
    totalMoves_++;
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
    
    // Generate translation
    Vector3 displacement = generateTranslationVector(2.0);  // Max 2 Angstrom move
    Vector3 newPos = oldPos + displacement;
    applyPeriodicBoundary(newPos);
    
    // Update position
    updateFragmentPosition(residueIdx, newPos);
    
    // Calculate energy after move
    result.energyAfter = calculateFragmentEnergy(residueIdx);
    result.deltaE = result.energyAfter - result.energyBefore;
    
    // Accept or reject
    bool accept = acceptMove(result.deltaE, 1.0, temperature_);
    
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
    
    // Generate rotation
    Quaternion rotation = generateRotationQuaternion(0.5);  // Max 0.5 radian rotation
    // Would multiply quaternions: newOrient = oldOrient * rotation
    Quaternion newOrient = oldOrient;
    // In a full implementation, would apply rotation here
    (void)rotation;  // Suppress unused variable warning
    newOrient.normalize();
    
    // Update orientation
    updateFragmentOrientation(residueIdx, newOrient);
    
    // Calculate energy after rotation
    result.energyAfter = calculateFragmentEnergy(residueIdx);
    result.deltaE = result.energyAfter - result.energyBefore;
    
    // Accept or reject
    bool accept = acceptMove(result.deltaE, 1.0, temperature_);
    
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
        return Vector3(uniform_(rng_) * 100 - 50,
                      uniform_(rng_) * 100 - 50,
                      uniform_(rng_) * 100 - 50);
    }
    
    return Vector3(uniform_(rng_) * state_->periodicBox[0] - state_->periodicBox[0]/2,
                  uniform_(rng_) * state_->periodicBox[1] - state_->periodicBox[1]/2,
                  uniform_(rng_) * state_->periodicBox[2] - state_->periodicBox[2]/2);
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
    if (energyCache_.valid) {
        return energyCache_.totalEnergy;
    }
    
    double totalEnergy = 0.0;
    
    if (!reservoir_) return totalEnergy;
    
    std::vector<int> instances = reservoir_->getActiveInstances();
    for (int idx : instances) {
        totalEnergy += calculateFragmentEnergy(idx);
    }
    
    // Divide by 2 to avoid double counting pairwise interactions
    totalEnergy /= 2.0;
    
    energyCache_.totalEnergy = totalEnergy;
    energyCache_.valid = true;
    
    return totalEnergy;
}

// Calculate fragment energy
double GCMCEngine::calculateFragmentEnergy(int residueIdx) {
    if (!reservoir_) return 0.0;
    
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (!instance) return 0.0;
    
    // Return cached energy if available
    if (instance->lastEnergyUpdate == totalMoves_) {
        return instance->energy_total;
    }
    
    double energy = calculateInteractionEnergy(residueIdx);
    
    // Update cache
    instance->energy_total = energy;
    instance->lastEnergyUpdate = totalMoves_;
    
    return energy;
}

// Calculate interaction energy
double GCMCEngine::calculateInteractionEnergy(int residueIdx) {
    double energy = 0.0;
    
    if (!reservoir_) return energy;
    
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (!instance) return energy;
    
    // Calculate interactions with all other fragments
    std::vector<int> others = reservoir_->getActiveInstances();
    for (int otherIdx : others) {
        if (otherIdx == residueIdx) continue;
        energy += calculatePairEnergy(residueIdx, otherIdx);
    }
    
    return energy;
}

// Calculate pair energy
double GCMCEngine::calculatePairEnergy(int idx1, int idx2) {
    if (!reservoir_) return 0.0;
    
    FragmentInstance* inst1 = reservoir_->getInstance(idx1);
    FragmentInstance* inst2 = reservoir_->getInstance(idx2);
    
    if (!inst1 || !inst2) return 0.0;
    
    double dist = minimumImageDistance(inst1->position, inst2->position);
    
    if (dist > cutoff_) return 0.0;
    
    // Simple LJ potential
    double sigma = 3.0;  // Angstrom
    double epsilon = 0.5;  // kJ/mol
    
    double r6 = std::pow(sigma / dist, 6);
    double r12 = r6 * r6;
    
    return 4.0 * epsilon * (r12 - r6);
}

// Calculate insertion bias
double GCMCEngine::calculateInsertionBias(const FragmentTemplate& tmpl,
                                         const Vector3& position,
                                         const Quaternion& orientation) {
    // Suppress unused parameter warnings
    (void)tmpl;
    (void)orientation;
    
    double bias = 1.0;
    
    if (cavityManager_) {
        // Convert to movement::Vector3
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
    if (!reservoir_) return 1.0;
    
    FragmentInstance* instance = reservoir_->getInstance(residueIdx);
    if (!instance) return 1.0;
    
    // Reverse of insertion bias
    double bias = 1.0;
    
    if (cavityManager_) {
        // Note: instance->position is already movement::Vector3
        double cavityScore = cavityManager_->getCavityScore(instance->position);
        if (cavityScore > 0) {
            bias /= cavityScore;
        }
    }
    
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

// Update fragment position
void GCMCEngine::updateFragmentPosition(int residueIdx, const Vector3& newPos) {
    if (!reservoir_) return;
    
    reservoir_->updatePosition(residueIdx, newPos);
    energyCache_.invalidate();
}

// Update fragment orientation
void GCMCEngine::updateFragmentOrientation(int residueIdx, const Quaternion& newOrient) {
    if (!reservoir_) return;
    
    reservoir_->updateOrientation(residueIdx, newOrient);
    energyCache_.invalidate();
}

// Apply periodic boundary conditions
void GCMCEngine::applyPeriodicBoundary(Vector3& position) {
    if (!state_) return;
    
    if (position.x < -state_->periodicBox[0]/2) position.x += state_->periodicBox[0];
    if (position.x > state_->periodicBox[0]/2) position.x -= state_->periodicBox[0];
    if (position.y < -state_->periodicBox[1]/2) position.y += state_->periodicBox[1];
    if (position.y > state_->periodicBox[1]/2) position.y -= state_->periodicBox[1];
    if (position.z < -state_->periodicBox[2]/2) position.z += state_->periodicBox[2];
    if (position.z > state_->periodicBox[2]/2) position.z -= state_->periodicBox[2];
}

// Calculate minimum image distance
double GCMCEngine::minimumImageDistance(const Vector3& r1, const Vector3& r2) {
    if (!state_) {
        return (r1 - r2).norm();
    }
    
    Vector3 dr = r1 - r2;
    
    // Apply minimum image convention
    if (dr.x > state_->periodicBox[0]/2) dr.x -= state_->periodicBox[0];
    if (dr.x < -state_->periodicBox[0]/2) dr.x += state_->periodicBox[0];
    if (dr.y > state_->periodicBox[1]/2) dr.y -= state_->periodicBox[1];
    if (dr.y < -state_->periodicBox[1]/2) dr.y += state_->periodicBox[1];
    if (dr.z > state_->periodicBox[2]/2) dr.z -= state_->periodicBox[2];
    if (dr.z < -state_->periodicBox[2]/2) dr.z += state_->periodicBox[2];
    
    return dr.norm();
}

} // namespace gcmc
} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc