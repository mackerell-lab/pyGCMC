// MultiInsertionCBMC.cpp
#include "../../energy/EnergyModule.hpp"
// Multi-insertion CBMC implementation

#include "MultiInsertionCBMC.hpp"
#include "../common/MovementUtils.hpp"
#include "../bias/CavityBias.hpp"
#include "../../energy/common/EnergyDirectCore.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <chrono>
#include <iostream>
#include <cstdlib>

#ifdef PYGCMC_USE_OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

using namespace model::montecarlo;
using gcmc::GCMCAcceptance;

MultiInsertionCBMC::MultiInsertionCBMC(const MultiInsertionConfig& config,
                                       CavityManager* cavityManager)
    : config_(config),
      rng_(static_cast<unsigned>(std::chrono::steady_clock::now().time_since_epoch().count())),
      uniform_(0.0, 1.0),
      cavityManager_(cavityManager) {
    
    // When using cavity bias, detailed balance is not guaranteed with Vbox normalization
    // Only show warning if PYGCMC_VERBOSE environment variable is set
    if (config_.useCavityBias) {
        static bool warningShown = false;
        if (!warningShown && std::getenv("PYGCMC_VERBOSE")) {
            std::cerr << "Note: Cavity bias enabled with Vbox normalization. "
                      << "Detailed balance with cavity-biased proposals is approximate. "
                      << "Use for enhanced sampling, not for strict DB validation." << std::endl;
            warningShown = true;
        }
    }
    
    resetStatistics();
}

void MultiInsertionCBMC::setSeed(uint64_t seed) {
    if (seed == 0) {
        rng_.seed(static_cast<unsigned>(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    } else {
        rng_.seed(static_cast<unsigned>(seed));
    }
}

std::pair<int, std::vector<InsertionRegion>> MultiInsertionCBMC::performMultiInsertion(
    MCState& state,
    int moleculeType,
    const MovementParams& params) {
    
    // Guard seeding: ensure RNG is seeded on first use
    if (stats_.totalAttempts == 0 && params.seed != 0) {
        setSeed(params.seed);
    }
    
    stats_.totalAttempts++;
    
    // 1. Divide box into regions
    auto allRegions = divideBoxIntoRegions(state);
    if (allRegions.empty()) {
        return {0, {}};
    }
    
    // 2. Select non-adjacent regions
    int numToSelect = std::min(config_.maxParallelInsertions, static_cast<int>(allRegions.size()));
    auto selectedRegions = selectNonAdjacentRegions(allRegions, numToSelect);
    
    // 3. Generate trial configurations
    generateTrialConfigurations(selectedRegions, moleculeType, state);
    
    // 4. Calculate energies
    batchCalculateEnergies(selectedRegions, state);
    
    // 5. Select optimal configurations
    selectOptimalConfigurations(selectedRegions, params.beta);
    
    // 6. Accept insertions
    auto acceptedRegions = acceptInsertions(
        selectedRegions,
        state,
        params.beta,
        params.chemicalPotential,
        params.thermalLambdaNm,
        moleculeType);
    
    // Count accepts
    int acceptCount = 0;
    for (const auto& region : acceptedRegions) {
        if (region.accepted) {
            acceptCount++;
        }
    }
    stats_.totalAccepts += acceptCount;
    
    return {acceptCount, acceptedRegions};
}

double MultiInsertionCBMC::getAcceptanceRate() const {
    if (stats_.totalAttempts == 0) return 0.0;
    return static_cast<double>(stats_.totalAccepts) / stats_.totalAttempts;
}

double MultiInsertionCBMC::getParallelEfficiency() const {
    return getAcceptanceRate();
}

void MultiInsertionCBMC::resetStatistics() {
    stats_ = Statistics();
}

std::vector<InsertionRegion> MultiInsertionCBMC::divideBoxIntoRegions(const MCState& state) {
    std::vector<InsertionRegion> regions;
    
    // Priority: use cavity points if available
    if (config_.useCavityBias && cavityManager_) {
        auto cavities = cavityManager_->findCavities(state);
        if (!cavities.empty()) {
            double minSepNm = config_.minSeparation; // Already in nm
            double regionRadius = 0.5 * minSepNm;     // Consistent with non-adjacent rules
            regions.reserve(cavities.size());
            
            for (const auto& pos : cavities) {
                InsertionRegion r;
                r.center[0] = pos.x;
                r.center[1] = pos.y;
                r.center[2] = pos.z;
                r.radius = regionRadius;
                r.gridIndex[0] = 0;  // Not used for cavity-based regions
                r.gridIndex[1] = 0;
                r.gridIndex[2] = 0;
                regions.push_back(r);
            }
            return regions;
        }
        // If no cavities found, fall back to uniform division
    }
    
    // Fallback: uniform grid division (original implementation)
    // Note: state.info.box is in nm
    double boxX = state.info.box[0];  // nm
    double boxY = state.info.box[1];  // nm
    double boxZ = state.info.box[2];  // nm
    
    double cellSize = config_.minSeparation;  // nm
    
    int nx = std::max(1, static_cast<int>(boxX / cellSize));
    int ny = std::max(1, static_cast<int>(boxY / cellSize));
    int nz = std::max(1, static_cast<int>(boxZ / cellSize));
    
    double actualCellX = boxX / nx;  // nm
    double actualCellY = boxY / ny;  // nm
    double actualCellZ = boxZ / nz;  // nm
    
    for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int iz = 0; iz < nz; ++iz) {
                InsertionRegion region;
                region.center[0] = (ix + 0.5) * actualCellX;
                region.center[1] = (iy + 0.5) * actualCellY;
                region.center[2] = (iz + 0.5) * actualCellZ;
                region.radius = std::min({actualCellX, actualCellY, actualCellZ}) * 0.5;
                region.gridIndex[0] = ix;
                region.gridIndex[1] = iy;
                region.gridIndex[2] = iz;
                regions.push_back(region);
            }
        }
    }
    
    return regions;
}

std::vector<InsertionRegion> MultiInsertionCBMC::selectNonAdjacentRegions(
    const std::vector<InsertionRegion>& allRegions,
    int numToSelect) {
    
    std::vector<InsertionRegion> selected;
    if (allRegions.empty() || numToSelect <= 0) {
        return selected;
    }
    
    // Create random order
    std::vector<size_t> indices(allRegions.size());
    for (size_t i = 0; i < allRegions.size(); ++i) {
        indices[i] = i;
    }
    std::shuffle(indices.begin(), indices.end(), rng_);
    
    // Calculate effective minimum separation for independence
    double minSepNm = config_.minSeparation;  // Already in nm
    
    if (config_.enforceIndependence) {
        // Use cutoff from config (should be set from MCState when initialized)
        double cutoffNm = config_.cutoffNm;
        
        // Calculate required separation for true independence
        // Need: separation > cutoff + sqrt(3)*disp + 2*molExtent
        double disp = (minSepNm * 0.5) * config_.displacementFraction;  // Upper bound estimate
        double requiredNm = cutoffNm + std::sqrt(3.0) * disp + 2.0 * config_.moleculeExtentNm;
        
        // Use the larger of user-specified and required separation
        double effectiveMinSep = std::max(minSepNm, requiredNm);
        
        // Log warning if user separation is too small
        if (minSepNm < requiredNm) {
            // Note: In production, should use proper logging
            // std::cerr << "Warning: minSeparation=" << minSepNm 
            //           << " nm is less than required=" << requiredNm 
            //           << " nm for independence\n";
        }
        
        minSepNm = effectiveMinSep;
    }
    
    for (size_t idx : indices) {
        const auto& candidate = allRegions[idx];
        
        // Check distance to all selected regions
        bool tooClose = false;
        for (const auto& sel : selected) {
            double dx = candidate.center[0] - sel.center[0];
            double dy = candidate.center[1] - sel.center[1];
            double dz = candidate.center[2] - sel.center[2];
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            
            if (dist < minSepNm) {
                tooClose = true;
                break;
            }
        }
        
        if (!tooClose) {
            selected.push_back(candidate);
            if (static_cast<int>(selected.size()) >= numToSelect) {
                break;
            }
        }
    }
    
    return selected;
}

void MultiInsertionCBMC::generateTrialConfigurations(
    std::vector<InsertionRegion>& regions,
    int moleculeType,
    const MCState& state) {
    
    for (auto& region : regions) {
        region.trialConfigs.clear();
        region.trialEnergies.clear();
        
        // Generate K trial configurations
        for (int k = 0; k < config_.numTrialsPerRegion; ++k) {
            // Create molecule at region center
            auto atoms = createMolecule(moleculeType, region.center);
            
            // Add random rotation and small displacement
            if (k > 0) {
                // Random rotation
                auto quat = generateRandomQuaternion();
                rotateWithQuaternion(atoms, quat);
                
                // Small random displacement within region
                double disp = region.radius * config_.displacementFraction;
                double dx = (uniform_(rng_) - 0.5) * disp;
                double dy = (uniform_(rng_) - 0.5) * disp;
                double dz = (uniform_(rng_) - 0.5) * disp;
                
                for (auto& atom : atoms) {
                    atom.x += dx;
                    atom.y += dy;
                    atom.z += dz;
                    
                    // Apply PBC
                    while (atom.x < 0) atom.x += state.info.box[0];
                    while (atom.x >= state.info.box[0]) atom.x -= state.info.box[0];
                    while (atom.y < 0) atom.y += state.info.box[1];
                    while (atom.y >= state.info.box[1]) atom.y -= state.info.box[1];
                    while (atom.z < 0) atom.z += state.info.box[2];
                    while (atom.z >= state.info.box[2]) atom.z -= state.info.box[2];
                }
            }
            
            region.trialConfigs.push_back(atoms);
        }
        
        // Initialize energies
        region.trialEnergies.resize(region.trialConfigs.size(), 0.0);
    }
}

// Helper: compute energies for a single region using a local state copy
static void calculateEnergiesForRegionImpl(InsertionRegion& region, MCState& state) {
    for (size_t t = 0; t < region.trialConfigs.size(); ++t) {
        const auto& trialAtoms = region.trialConfigs[t];
        int startIdx = state.activeAtomCount;
        for (const auto& atom : trialAtoms) {
            state.addAtom(atom);
        }
        
        MCResidue tempRes;
        tempRes.atomStart = startIdx;
        tempRes.atomCount = static_cast<int>(trialAtoms.size());
        tempRes.type = 0;
        tempRes.active = true;
        int resIdx = state.addResidue(tempRes);
        
        // Compute only the interaction energy for the new residue
        platform::cpu::computeResidueEnergyCutoffPBC(state, resIdx);
        region.trialEnergies[t] = state.residues[resIdx].energy_vdw + 
                                  state.residues[resIdx].energy_elec;
        
        state.removeResidue(resIdx);
        for (size_t i = 0; i < trialAtoms.size(); ++i) {
            state.removeAtom(state.activeAtomCount - 1);
        }
    }
}

void MultiInsertionCBMC::batchCalculateEnergies(
    std::vector<InsertionRegion>& regions,
    MCState& state) {
    
#ifdef PYGCMC_USE_OPENMP
    // Use static schedule for deterministic ordering
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < static_cast<int>(regions.size()); ++i) {
        MCState localState = state; // Thread-local copy
        calculateEnergiesForRegionImpl(regions[i], localState);
    }
#else
    // Serial execution
    for (auto& region : regions) {
        calculateEnergiesForRegionImpl(region, state);
    }
#endif
}

void MultiInsertionCBMC::selectOptimalConfigurations(
    std::vector<InsertionRegion>& regions,
    double beta) {
    
    for (auto& region : regions) {
        // Calculate Rosenbluth weight
        region.rosenbluthWeight = calculateRosenbluthWeight(region.trialEnergies, beta);
        
        // Select configuration using Gumbel-max trick
        std::vector<double> logWeights;
        for (double energy : region.trialEnergies) {
            logWeights.push_back(-beta * energy);
        }
        
        region.selectedConfig = selectByGumbelMax(logWeights);
    }
}

std::vector<InsertionRegion> MultiInsertionCBMC::acceptInsertions(
    std::vector<InsertionRegion>& regions,
    const MCState& state,
    double beta,
    double chemicalPotential,
    double lambdaNm,
    int moleculeType) {
    
    // Count current molecules
    int currentN = state.activeResidueCount;
    
    // Precompute box volume in nm^3
    double Vbox = state.info.box[0] * state.info.box[1] * state.info.box[2];
    
    double lambda = (lambdaNm > 0.0) ? lambdaNm : 1.0;
    double logLambda3 = 3.0 * std::log(lambda);

    // Process each region independently
    for (auto& region : regions) {
        // Find max log weight for numerical stability
        double maxLogW = -std::numeric_limits<double>::infinity();
        int Keff = 0;
        
        for (double e : region.trialEnergies) {
            double lw = -beta * e;
            if (!std::isfinite(lw)) continue;
            if (lw > maxLogW) maxLogW = lw;
        }
        
        // If all energies are infinite, reject
        if (!std::isfinite(maxLogW)) {
            region.accepted = false;
            region.acceptanceProbability = 0.0;
            continue;
        }
        
        // Compute Rosenbluth weight
        double sumExp = 0.0;
        
        for (double e : region.trialEnergies) {
            double lw = -beta * e;
            if (!std::isfinite(lw)) continue;
            sumExp += std::exp(lw - maxLogW);
            Keff++;
        }
        
        if (Keff == 0) {
            region.accepted = false;
            region.acceptanceProbability = 0.0;
            continue;
        }
        
        double logW = maxLogW + std::log(sumExp);
        
        region.countBefore = currentN;
        region.effectiveTrials = Keff;
        if (region.selectedConfig >= 0 &&
            region.selectedConfig < static_cast<int>(region.trialEnergies.size())) {
            region.selectedEnergy = region.trialEnergies[region.selectedConfig];
        } else {
            region.selectedEnergy = 0.0;
        }

        double logVolume = std::log(std::max(1e-30, Vbox));
        double logKeff = std::log(static_cast<double>(std::max(Keff, 1)));
        double logRosen = logW - logKeff;
        region.rosenbluthWeight = std::exp(logRosen);

        region.grandTerms.speciesId = moleculeType;
        region.grandTerms.countBefore = region.countBefore;
        region.grandTerms.countAfter = region.countBefore + 1;
        region.grandTerms.beta = beta;
        region.grandTerms.chemicalPotential = chemicalPotential;
        region.grandTerms.deltaEnergy = region.selectedEnergy;
        region.grandTerms.logVolume = logVolume;
        region.grandTerms.logLambda3 = logLambda3;
        region.grandTerms.logProposalForward = 0.0;
        region.grandTerms.logProposalReverse = 0.0;
        region.grandTerms.logCavityForward = 0.0;
        region.grandTerms.logCavityReverse = 0.0;
        region.grandTerms.logRosenbluthForward = logRosen;
        region.grandTerms.logRosenbluthReverse = 0.0;
        region.grandTerms.logExtraForward = 0.0;
        region.grandTerms.logExtraReverse = 0.0;

        region.grandEval = GCMCAcceptance::evaluate(
            region.grandTerms,
            GCMCAcceptance::MoveType::INSERTION);
        region.logAcceptanceRatio = region.grandEval.logRatio;
        region.acceptanceProbability = region.grandEval.probability;

        // Accept or reject
        region.accepted = (uniform_(rng_) < region.acceptanceProbability);
        
        // Update molecule count if accepted
        if (region.accepted) {
            currentN++;
        }
    }
    
    return regions;
}

std::vector<MCAtom> MultiInsertionCBMC::createMolecule(
    int moleculeType,
    const std::array<double, 3>& position) {
    
    (void)moleculeType;  // Suppress unused parameter warning
    std::vector<MCAtom> atoms;
    
    // Create water molecule (TIP3P)
    // Note: sigma and epsilon are stored in force field parameters, not in MCAtom
    MCAtom oxygen;
    oxygen.x = position[0];
    oxygen.y = position[1];
    oxygen.z = position[2];
    oxygen.charge = -0.834f;
    oxygen.type = 0;  // Oxygen type
    atoms.push_back(oxygen);
    
    MCAtom hydrogen1;
    hydrogen1.x = position[0] + 0.0957f;
    hydrogen1.y = position[1];
    hydrogen1.z = position[2];
    hydrogen1.charge = 0.417f;
    hydrogen1.type = 1;  // Hydrogen type
    atoms.push_back(hydrogen1);
    
    MCAtom hydrogen2;
    hydrogen2.x = position[0] - 0.0239f;
    hydrogen2.y = position[1] + 0.0927f;
    hydrogen2.z = position[2];
    hydrogen2.charge = 0.417f;
    hydrogen2.type = 1;  // Hydrogen type
    atoms.push_back(hydrogen2);
    
    return atoms;
}

std::array<double, 4> MultiInsertionCBMC::generateRandomQuaternion() {
    // Generate random quaternion for uniform rotation
    double u1 = uniform_(rng_);
    double u2 = uniform_(rng_);
    double u3 = uniform_(rng_);
    
    double q0 = std::sqrt(1 - u1) * std::sin(2 * M_PI * u2);
    double q1 = std::sqrt(1 - u1) * std::cos(2 * M_PI * u2);
    double q2 = std::sqrt(u1) * std::sin(2 * M_PI * u3);
    double q3 = std::sqrt(u1) * std::cos(2 * M_PI * u3);
    
    return {q0, q1, q2, q3};
}

void MultiInsertionCBMC::rotateWithQuaternion(
    std::vector<MCAtom>& atoms,
    const std::array<double, 4>& quat) {
    
    if (atoms.empty()) return;
    
    // Calculate center of mass
    double cx = 0, cy = 0, cz = 0;
    for (const auto& atom : atoms) {
        cx += atom.x;
        cy += atom.y;
        cz += atom.z;
    }
    cx /= atoms.size();
    cy /= atoms.size();
    cz /= atoms.size();
    
    // Rotate around center
    double q0 = quat[0], q1 = quat[1], q2 = quat[2], q3 = quat[3];
    
    for (auto& atom : atoms) {
        // Translate to origin
        double x = atom.x - cx;
        double y = atom.y - cy;
        double z = atom.z - cz;
        
        // Apply rotation
        double xx = (q0*q0 + q1*q1 - q2*q2 - q3*q3) * x + 2*(q1*q2 - q0*q3) * y + 2*(q1*q3 + q0*q2) * z;
        double yy = 2*(q1*q2 + q0*q3) * x + (q0*q0 - q1*q1 + q2*q2 - q3*q3) * y + 2*(q2*q3 - q0*q1) * z;
        double zz = 2*(q1*q3 - q0*q2) * x + 2*(q2*q3 + q0*q1) * y + (q0*q0 - q1*q1 - q2*q2 + q3*q3) * z;
        
        // Translate back
        atom.x = xx + cx;
        atom.y = yy + cy;
        atom.z = zz + cz;
    }
}

double MultiInsertionCBMC::calculateRosenbluthWeight(
    const std::vector<double>& energies,
    double beta) {
    
    double weight = 0.0;
    for (double energy : energies) {
        weight += std::exp(-beta * energy);
    }
    return weight;
}

int MultiInsertionCBMC::selectByGumbelMax(const std::vector<double>& logWeights) {
    if (logWeights.empty()) return -1;
    
    int bestIdx = 0;
    double bestScore = -std::numeric_limits<double>::infinity();
    
    for (size_t i = 0; i < logWeights.size(); ++i) {
        // Add Gumbel noise
        double u = uniform_(rng_);
        double gumbel = -std::log(-std::log(u + 1e-10) + 1e-10);
        double score = logWeights[i] + gumbel;
        
        if (score > bestScore) {
            bestScore = score;
            bestIdx = static_cast<int>(i);
        }
    }
    
    return bestIdx;
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc
