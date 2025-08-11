#include "DrudeSCFOM.hpp"
#include "../DrudeStructures.hpp"  // Contains DrudeConstants
#include <cmath>
#include <algorithm>
#include <iostream>
#include <map>
#include <unordered_set>
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace exp {

void DrudeSCFOM::applyPBC(double& dx, double& dy, double& dz, 
                               const std::array<double, 3>& box) {
    if (box[0] > 0) dx -= box[0] * std::round(dx / box[0]);
    if (box[1] > 0) dy -= box[1] * std::round(dy / box[1]);
    if (box[2] > 0) dz -= box[2] * std::round(dz / box[2]);
}

double DrudeSCFOM::tholeS1(double r, double alpha_i, double alpha_j, double thole_pair) const {
    // CHARMM/OpenMM standard S1 screening function
    // S1(u) = 1 - (1 + u/2)exp(-u)
    // where u = a*r / (alpha_i * alpha_j)^(1/6)
    
    // If thole parameter is zero or negligible, return 1.0 (no screening)
    if (thole_pair <= 1e-10 || alpha_i <= 1e-14 || alpha_j <= 1e-14) {
        return 1.0;  // No screening
    }
    
    double alpha_eff = std::pow(alpha_i * alpha_j, 1.0/6.0);
    double u = thole_pair * r / alpha_eff;
    
    if (u > 50.0) {
        return 1.0;  // Screening is negligible for large u
    }
    
    double exp_u = std::exp(-u);
    // S1 function for point charge screening (CHARMM standard)
    return 1.0 - (1.0 + 0.5 * u) * exp_u;
}

/* DEPRECATED: Non-standard S3 screening function
 * The S3 function is not part of the standard CHARMM/OpenMM Drude model.
 * Use S1 screening with point charges instead.
 */
/*
double DrudeSCFOM::tholeS3(double r, double alpha_i, double alpha_j, double thole_pair) const {
    // Thole S3 screening for 1/r^3 dipole-dipole interactions
    // S3(u) = 1 - exp(-u) * (1 + u + u^2/2)
    // where u = a*r / (alpha_i * alpha_j)^(1/6)
    // and a is the pair-level Thole parameter
    
    // If thole parameter is zero or negligible, return 1.0 (no screening)
    if (thole_pair <= 1e-10 || alpha_i <= 1e-14 || alpha_j <= 1e-14) {
        return 1.0;  // No screening
    }
    
    double alpha_eff = std::pow(alpha_i * alpha_j, 1.0/6.0);
    double u = thole_pair * r / alpha_eff;
    
    if (u > 50.0) {
        return 1.0;  // Screening is negligible for large u
    }
    
    double exp_u = std::exp(-u);
    // S3 function for 1/r^3 dipole field screening
    return 1.0 - exp_u * (1.0 + u + 0.5 * u * u);
}
*/

/* DEPRECATED: Non-standard S5 screening function  
 * The S5 function is not part of the standard CHARMM/OpenMM Drude model.
 * Use S1 screening with point charges instead.
 */
/*
double DrudeSCFOM::tholeS5(double r, double alpha_i, double alpha_j, double thole_pair) const {
    // Thole S5 screening for 1/r^5 tensor component
    // S5(u) = 1 - exp(-u) * (1 + u + u^2/2 + u^3/6)
    // where u = a*r / (alpha_i * alpha_j)^(1/6)
    // and a is the pair-level Thole parameter
    
    // If thole parameter is zero or negligible, return 1.0 (no screening)
    if (thole_pair <= 1e-10 || alpha_i <= 1e-14 || alpha_j <= 1e-14) {
        return 1.0;  // No screening
    }
    
    double alpha_eff = std::pow(alpha_i * alpha_j, 1.0/6.0);
    double u = thole_pair * r / alpha_eff;
    
    if (u > 50.0) {
        return 1.0;  // Screening is negligible for large u
    }
    
    double exp_u = std::exp(-u);
    // S5 function for 1/r^5 tensor component screening
    return 1.0 - exp_u * (1.0 + u + 0.5 * u * u + u * u * u / 6.0);
}
*/

bool DrudeSCFOM::optimize(model::MCState& state,
                              const std::vector<DrudeParticle>& particles,
                              const std::vector<ScreenedPair>& pairs,
                              const DrudeSCFParams& params) {
    if (particles.empty()) {
        return true;
    }
    
    // Initialize electric field vectors
    std::vector<Vec3> electricField(particles.size(), {0.0, 0.0, 0.0});
    
    // Best-so-far tracking
    double bestEnergy = 1e300;
    std::vector<DrudePosition> bestPositions(particles.size());
    
    // Adaptive damping
    double damping = params.dampingFactor;
    double previousEnergy = 1e300;
    int stallCount = 0;
    
    // DIIS acceleration data
    DIISData diis;
    diis.maxHistory = 5;
    diis.enabled = false;  // Enable after initial iterations
    
    // SCF iteration
    for (int iter = 0; iter < params.maxIterations; ++iter) {
        // Clear electric field
        std::fill(electricField.begin(), electricField.end(), Vec3{0.0, 0.0, 0.0});
        
        // Calculate total electric field
        calculateElectricField(state, particles, pairs, electricField, params);
        
        // Debug output for first few iterations
        if (iter < 5 || iter == params.maxIterations - 1) {
            double totalFieldMag = 0.0;
            for (const auto& field : electricField) {
                totalFieldMag += std::sqrt(field[0]*field[0] + field[1]*field[1] + field[2]*field[2]);
            }
            std::cout << "SCF iter " << iter << ": avg field magnitude = " << totalFieldMag/electricField.size() << std::endl;
        }
        
        // Check convergence based on forces (skip first few iterations to ensure proper convergence)
        if (iter > 2) {
            double maxForce = 0.0;
            for (size_t i = 0; i < particles.size(); ++i) {
                const auto& p = particles[i];
                if (p.polarizability < 1e-14) continue;  // Skip frozen particles
                
                const auto& drude = state.atoms[p.drudeIndex];
                const auto& parent = state.atoms[p.parentIndex];
                
                double dx = drude.x - parent.x;
                double dy = drude.y - parent.y;
                double dz = drude.z - parent.z;
                
                std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
                applyPBC(dx, dy, dz, box);
                
                // Force = q*E - k*d
                double fx = p.charge * electricField[i][0] - p.kSpring * dx;
                double fy = p.charge * electricField[i][1] - p.kSpring * dy;
                double fz = p.charge * electricField[i][2] - p.kSpring * dz;
                
                double force2 = fx*fx + fy*fy + fz*fz;
                maxForce = std::max(maxForce, std::sqrt(force2));
            }
            
            // Debug output
            if (iter < 5) {
                std::cout << "  iter " << iter << ": maxForce = " << maxForce 
                          << " (tolerance = " << params.tolerance << ")" << std::endl;
            }
            
            // Check convergence
            if (maxForce < params.tolerance) {
                std::cout << "SCF converged at iteration " << iter << " with maxForce = " << maxForce << std::endl;
                return true;  // Converged
            }
            
            // Debug warning if approaching max iterations
            if (iter == params.maxIterations - 1) {
                std::cout << "WARNING: SCF did not converge! Final maxForce = " << maxForce 
                          << " (tolerance = " << params.tolerance << ")" << std::endl;
            }
        }
        
        // Enable DIIS after initial iterations for stability
        if (iter == 5) {
            diis.enabled = true;
        }
        
        // Try DIIS acceleration if enabled
        bool diisApplied = false;
        if (diis.enabled) {
            diisApplied = applyDIIS(diis, particles, electricField, state);
        }
        
        // If DIIS not applied, use regular damped update
        if (!diisApplied) {
            double maxStep = 0.1;  // 0.1 nm max step per iteration (allow larger steps)
            double hardWall = params.enableHardWall ? params.maxDrudeDistance : 1e9;
            updateDrudePositions(state, particles, electricField, damping, maxStep, hardWall);
        }
        
        // Calculate current energy
        double currentEnergy = calculateSpringEnergy(state, particles);
        
        // Track best-so-far
        if (currentEnergy < bestEnergy) {
            bestEnergy = currentEnergy;
            for (size_t i = 0; i < particles.size(); ++i) {
                const auto& p = particles[i];
                const auto& drude = state.atoms[p.drudeIndex];
                const auto& parent = state.atoms[p.parentIndex];
                bestPositions[i].dx = drude.x - parent.x;
                bestPositions[i].dy = drude.y - parent.y;
                bestPositions[i].dz = drude.z - parent.z;
            }
            stallCount = 0;
        } else {
            stallCount++;
        }
        
        // Adaptive damping
        if (currentEnergy < previousEnergy) {
            damping *= 0.95;  // Reduce damping if improving
            damping = std::max(damping, 0.05);
        } else if (stallCount > 5) {
            damping *= 1.1;  // Increase damping if stalled
            damping = std::min(damping, 0.5);
            stallCount = 0;
        }
        
        previousEnergy = currentEnergy;
    }
    
    // Only revert to best-so-far if we had instability
    // Otherwise keep the last positions (which should be closer to self-consistency)
    if (bestEnergy < 1e100) {  // We found at least one good configuration
        // Check if current configuration is reasonable
        double finalEnergy = calculateSpringEnergy(state, particles);
        if (finalEnergy > bestEnergy * 2.0) {  // Current is much worse than best
            // Revert to best-so-far
            for (size_t i = 0; i < particles.size(); ++i) {
                const auto& p = particles[i];
                auto& drude = state.atoms[p.drudeIndex];
                const auto& parent = state.atoms[p.parentIndex];
                drude.x = parent.x + bestPositions[i].dx;
                drude.y = parent.y + bestPositions[i].dy;
                drude.z = parent.z + bestPositions[i].dz;
            }
        }
        // Otherwise keep current positions
    }
    
    return false;  // Did not converge
}

void DrudeSCFOM::calculateElectricField(const model::MCState& state,
                                            const std::vector<DrudeParticle>& particles,
                                            const std::vector<ScreenedPair>& pairs,
                                            std::vector<Vec3>& electricField,
                                            const DrudeSCFParams& params) const {
    // Field breakdown debugging: Track external and induced fields separately
    std::vector<Vec3> E_ext(particles.size(), {0.0, 0.0, 0.0});
    std::vector<Vec3> E_ind(particles.size(), {0.0, 0.0, 0.0});
    
    // Calculate external field (from fixed charges)
    calculateExternalField(state, particles, pairs, E_ext, params);
    
    // Calculate induced field (from dipoles with Thole screening)
    calculateInducedField(state, particles, pairs, E_ind);
    
    // Combine fields and debug output
    for (size_t i = 0; i < particles.size(); ++i) {
        electricField[i][0] = E_ext[i][0] + E_ind[i][0];
        electricField[i][1] = E_ext[i][1] + E_ind[i][1];
        electricField[i][2] = E_ext[i][2] + E_ind[i][2];
        
        // Debug: Print field breakdown for first particle
        if (i == 0 && !pairs.empty()) {
            double ext_mag = std::sqrt(E_ext[i][0]*E_ext[i][0] + E_ext[i][1]*E_ext[i][1] + E_ext[i][2]*E_ext[i][2]);
            double ind_mag = std::sqrt(E_ind[i][0]*E_ind[i][0] + E_ind[i][1]*E_ind[i][1] + E_ind[i][2]*E_ind[i][2]);
            double tot_mag = std::sqrt(electricField[i][0]*electricField[i][0] + 
                                      electricField[i][1]*electricField[i][1] + 
                                      electricField[i][2]*electricField[i][2]);
            std::cout << "FIELD_BREAKDOWN[0]: E_ext=" << ext_mag 
                      << " E_ind=" << ind_mag 
                      << " E_total=" << tot_mag 
                      << " (ind/ext=" << (ext_mag > 1e-10 ? ind_mag/ext_mag : 0) << ")" 
                      << std::endl;
        }
    }
}

void DrudeSCFOM::calculateExternalField(const model::MCState& state,
                                            const std::vector<DrudeParticle>& particles,
                                            const std::vector<ScreenedPair>& pairs,
                                            std::vector<Vec3>& electricField,
                                            const DrudeSCFParams& params) const {
    // Build mask for ALL Drude atoms to exclude them from external field
    // This prevents double-counting of induced-induced interactions
    std::vector<bool> isDrude(state.activeAtomCount, false);
    for (const auto& p : particles) {
        if (p.drudeIndex >= 0 && p.drudeIndex < state.activeAtomCount) {
            isDrude[p.drudeIndex] = true;
        }
    }
    
    // Build exclusion list if excludePartnerParentInExternalField is enabled
    std::vector<std::vector<int>> excludeParents(particles.size());
    if (params.excludePartnerParentInExternalField) {
        for (const auto& pair : pairs) {
            if (static_cast<size_t>(pair.dipole1) < particles.size() && 
                static_cast<size_t>(pair.dipole2) < particles.size()) {
                // For each Drude, exclude the parent of its partner to avoid double counting
                excludeParents[pair.dipole1].push_back(particles[pair.dipole2].parentIndex);
                excludeParents[pair.dipole2].push_back(particles[pair.dipole1].parentIndex);
            }
        }
    }
    
    // Calculate external field from non-Drude charges ONLY
    // IMPORTANT: Field should be calculated at parent position, not Drude position!
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        const auto& parentAtom = state.atoms[particle.parentIndex];  // Use parent position!
        
        Vec3 field = {0.0, 0.0, 0.0};
        
        for (int j = 0; j < state.activeAtomCount; ++j) {
            // Skip self and parent
            if (j == particle.drudeIndex || j == particle.parentIndex) continue;
            
            // CRITICAL: Skip ALL Drude charges to prevent double-counting
            if (isDrude[j]) continue;
            
            // CRITICAL: Skip all atoms in the same molecule (intramolecular exclusion)
            // This prevents self-polarization in water and other molecules
            if (inSameMolecule(particle.parentIndex, j, state)) {
                continue;
            }
            
            // If excludePartnerParentInExternalField is enabled, skip partner parents
            if (params.excludePartnerParentInExternalField) {
                if (std::find(excludeParents[i].begin(), excludeParents[i].end(), j) 
                    != excludeParents[i].end()) {
                    continue;
                }
            }
            
            const auto& atom = state.atoms[j];
            if (std::abs(atom.charge) < 1e-10) continue;
            
            // Calculate distance vector FROM charge TO parent (field direction)
            // Field points from source (atom) to target (parent)
            double dx = parentAtom.x - atom.x;
            double dy = parentAtom.y - atom.y;
            double dz = parentAtom.z - atom.z;
            
            std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
            applyPBC(dx, dy, dz, box);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1e-12) continue;  // Skip if too close
            
            double r = std::sqrt(r2);
            double r3 = r2 * r;
            
            // Electric field: E = k * q / r^2 * r_hat
            double factor = DrudeConstants::ONE_4PI_EPS0 * atom.charge / r3;
            
            field[0] += factor * dx;
            field[1] += factor * dy;
            field[2] += factor * dz;
        }
        
        electricField[i][0] += field[0];
        electricField[i][1] += field[1];
        electricField[i][2] += field[2];
        
        // Debug first particle external field
        if (i == 0) {
            std::cout << "EXT_FIELD[0]: E_x=" << field[0] << " E_y=" << field[1] 
                      << " E_z=" << field[2] << std::endl;
        }
    }
}

void DrudeSCFOM::calculateInducedField(const model::MCState& state,
                                           const std::vector<DrudeParticle>& particles,
                                           const std::vector<ScreenedPair>& pairs,
                                           std::vector<Vec3>& electricField) const {
    // Dispatch to appropriate algorithm
    switch (algorithm) {
        case DrudeAlgorithm::S1_POINT_CHARGE:
            calculateInducedFieldS1(state, particles, pairs, electricField);
            break;
        /* DEPRECATED: Non-standard S3/S5 model
        case DrudeAlgorithm::S3S5_DIPOLE_TENSOR:
            calculateInducedFieldS3S5(state, particles, pairs, electricField);
            break;
        case DrudeAlgorithm::S1_DIPOLE_FIELD:
            // TODO: Implement hybrid approach
            calculateInducedFieldS1(state, particles, pairs, electricField);
            break;
        */
        case DrudeAlgorithm::DIRECT_COULOMB:
            // No induced field for direct Coulomb (testing only)
            break;
        default:
            // Default to standard S1 model
            calculateInducedFieldS1(state, particles, pairs, electricField);
            break;
    }
}

void DrudeSCFOM::calculateInducedFieldS1(const model::MCState& state,
                                          const std::vector<DrudeParticle>& particles,
                                          const std::vector<ScreenedPair>& pairs,
                                          std::vector<Vec3>& electricField) const {
    // CHARMM/OpenMM standard: 4-point charge model with S1 screening
    // Each dipole is represented by parent(+q) and Drude(-q) charges
    // All 4 charge-charge interactions are computed with S1 screening
    
    // Build map of screened pairs for quick lookup
    std::map<std::pair<int,int>, double> screenedPairs;
    for (const auto& pair : pairs) {
        screenedPairs[{pair.dipole1, pair.dipole2}] = pair.thole;
        screenedPairs[{pair.dipole2, pair.dipole1}] = pair.thole;
    }
    
    // Process all dipole pairs
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            const auto& p1 = particles[i];
            const auto& p2 = particles[j];
            
            // Get Thole parameter (0 if not a screened pair)
            double a_pair = 0.0;
            auto it = screenedPairs.find({static_cast<int>(i), static_cast<int>(j)});
            if (it != screenedPairs.end()) {
                a_pair = it->second;
            }
            // ✅ Key fix: Process ALL dipole pairs, even when a_pair=0 (S(u)=1)
            
            const auto& parent1 = state.atoms[p1.parentIndex];
            const auto& drude1 = state.atoms[p1.drudeIndex];
            const auto& parent2 = state.atoms[p2.parentIndex];
            const auto& drude2 = state.atoms[p2.drudeIndex];
            
            std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
            
            // Calculate parent-parent distance for u calculation (OpenMM standard)
            double dx_pp = parent2.x - parent1.x;
            double dy_pp = parent2.y - parent1.y;
            double dz_pp = parent2.z - parent1.z;
            applyPBC(dx_pp, dy_pp, dz_pp, box);
            double r_pp = std::sqrt(dx_pp*dx_pp + dy_pp*dy_pp + dz_pp*dz_pp);
            
            // Four charge-charge interactions with S1 screening:
            // 1. Parent1 - Parent2
            // 2. Parent1 - Drude2
            // 3. Drude1 - Parent2
            // 4. Drude1 - Drude2
            
            // Helper lambda for screened Coulomb field calculation
            auto computeScreenedField = [&](const auto& atom1, const auto& atom2, 
                                           double q2, Vec3& field) {
                double dx = atom1.x - atom2.x;
                double dy = atom1.y - atom2.y;
                double dz = atom1.z - atom2.z;
                applyPBC(dx, dy, dz, box);
                
                double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 < 1e-12) return;
                
                double r = std::sqrt(r2);
                // Use parent-parent distance for u calculation (OpenMM standard)
                double s1 = tholeS1(r_pp, p1.polarizability, p2.polarizability, a_pair);
                
                // Field = k * q * S(r) / r^3 * r_vec
                double factor = DrudeConstants::ONE_4PI_EPS0 * q2 * s1 / (r2 * r);
                field[0] += factor * dx;
                field[1] += factor * dy;
                field[2] += factor * dz;
            };
            
            // Calculate polarization charges
            double qPolP1 = -p1.charge;  // Parent gets +|qD|
            double qPolD1 = p1.charge;    // Drude keeps -qD
            double qPolP2 = -p2.charge;
            double qPolD2 = p2.charge;
            
            // Field on dipole 1
            Vec3 field1 = {0, 0, 0};
            computeScreenedField(parent1, parent2, qPolP2, field1);  // P1 from P2
            computeScreenedField(parent1, drude2, qPolD2, field1);   // P1 from D2
            computeScreenedField(drude1, parent2, qPolP2, field1);   // D1 from P2
            computeScreenedField(drude1, drude2, qPolD2, field1);    // D1 from D2
            
            // Field on dipole 2
            Vec3 field2 = {0, 0, 0};
            computeScreenedField(parent2, parent1, qPolP1, field2);  // P2 from P1
            computeScreenedField(parent2, drude1, qPolD1, field2);   // P2 from D1
            computeScreenedField(drude2, parent1, qPolP1, field2);   // D2 from P1
            computeScreenedField(drude2, drude1, qPolD1, field2);    // D2 from D1
            
            // Add to total field (field acts on Drude particles)
            electricField[i][0] += field1[0];
            electricField[i][1] += field1[1];
            electricField[i][2] += field1[2];
            
            electricField[j][0] += field2[0];
            electricField[j][1] += field2[1];
            electricField[j][2] += field2[2];
        }
    }
}

/* DEPRECATED: Non-standard S3/S5 dipole tensor model
 * This implementation uses a non-standard dipole tensor field approach
 * that is not consistent with CHARMM/OpenMM standards.
 * Use S1_POINT_CHARGE algorithm instead.
 */
/*
void DrudeSCFOM::calculateInducedFieldS3S5(const model::MCState& state,
                                           const std::vector<DrudeParticle>& particles,
                                           const std::vector<ScreenedPair>& pairs,
                                           std::vector<Vec3>& electricField) const {
    // Calculate induced field using dipole tensor model with S3/S5 Thole screening
    // Field: E = k/r³ [3 S5(u) (μ·n) n − S3(u) μ]
    // where μ = q_D * (D - P) is the dipole moment based on current Drude positions
    
    // Build map of screened pairs for quick lookup
    std::map<std::pair<int,int>, double> screenedPairs;
    for (const auto& pair : pairs) {
        screenedPairs[{pair.dipole1, pair.dipole2}] = pair.thole;
        screenedPairs[{pair.dipole2, pair.dipole1}] = pair.thole;
    }
    
    // Process ALL dipole pairs
    for (size_t i = 0; i < particles.size(); ++i) {
        for (size_t j = i + 1; j < particles.size(); ++j) {
            const auto& p1 = particles[i];
            const auto& p2 = particles[j];
            const auto& parent1 = state.atoms[p1.parentIndex];
            const auto& drude1 = state.atoms[p1.drudeIndex];
            const auto& parent2 = state.atoms[p2.parentIndex];
            const auto& drude2 = state.atoms[p2.drudeIndex];
            
            std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
            
            // Calculate current dipole moments μ = q_D * (D - P)
            double dx1 = drude1.x - parent1.x;
            double dy1 = drude1.y - parent1.y;
            double dz1 = drude1.z - parent1.z;
            applyPBC(dx1, dy1, dz1, box);
            
            double dx2 = drude2.x - parent2.x;
            double dy2 = drude2.y - parent2.y;
            double dz2 = drude2.z - parent2.z;
            applyPBC(dx2, dy2, dz2, box);
            
            Vec3 mu1 = {p1.charge * dx1, p1.charge * dy1, p1.charge * dz1};
            Vec3 mu2 = {p2.charge * dx2, p2.charge * dy2, p2.charge * dz2};
            
            // Skip if both dipoles are zero (no induced field)
            double mu1_mag2 = mu1[0]*mu1[0] + mu1[1]*mu1[1] + mu1[2]*mu1[2];
            double mu2_mag2 = mu2[0]*mu2[0] + mu2[1]*mu2[1] + mu2[2]*mu2[2];
            if (mu1_mag2 < 1e-20 && mu2_mag2 < 1e-20) continue;
            
            // Calculate distance between parent atoms (field calculation points)
            double rx = parent2.x - parent1.x;
            double ry = parent2.y - parent1.y;
            double rz = parent2.z - parent1.z;
            applyPBC(rx, ry, rz, box);
            
            double r2 = rx*rx + ry*ry + rz*rz;
            if (r2 < 1e-12) continue;
            
            double r = std::sqrt(r2);
            double invr3 = 1.0 / (r2 * r);
            
            // Unit vector from dipole 1 to dipole 2
            double nx = rx / r;
            double ny = ry / r;
            double nz = rz / r;
            
            // Check if this pair has Thole screening
            double a_pair = 0.0;
            auto it = screenedPairs.find({static_cast<int>(i), static_cast<int>(j)});
            if (it != screenedPairs.end()) {
                a_pair = it->second;  // Use pair.thole directly, no multiplication
            }
            
            // Calculate S3 and S5 screening functions
            double s3 = tholeS3(r, p1.polarizability, p2.polarizability, a_pair);
            double s5 = tholeS5(r, p1.polarizability, p2.polarizability, a_pair);
            
            // Debug output for first pair
            static bool debugPrinted = false;
            if (!debugPrinted && a_pair > 0) {
                std::cout << "THOLE_DEBUG: r=" << r << " alpha1=" << p1.polarizability 
                          << " alpha2=" << p2.polarizability << " a_pair=" << a_pair 
                          << " s3=" << s3 << " s5=" << s5 << std::endl;
                debugPrinted = true;
            }
            
            // Prefactor
            double prefactor = DrudeConstants::ONE_4PI_EPS0 * invr3;
            
            // Field on dipole 1 from dipole 2: E1 = k/r³ [3 S5(u) (μ2·n) n − S3(u) μ2]
            double mu2_dot_n = mu2[0]*nx + mu2[1]*ny + mu2[2]*nz;
            Vec3 E1 = {
                prefactor * (3.0 * s5 * mu2_dot_n * nx - s3 * mu2[0]),
                prefactor * (3.0 * s5 * mu2_dot_n * ny - s3 * mu2[1]),
                prefactor * (3.0 * s5 * mu2_dot_n * nz - s3 * mu2[2])
            };
            
            // Field on dipole 2 from dipole 1: symmetric formula
            double mu1_dot_n = mu1[0]*nx + mu1[1]*ny + mu1[2]*nz;
            Vec3 E2 = {
                prefactor * (3.0 * s5 * mu1_dot_n * nx - s3 * mu1[0]),
                prefactor * (3.0 * s5 * mu1_dot_n * ny - s3 * mu1[1]),
                prefactor * (3.0 * s5 * mu1_dot_n * nz - s3 * mu1[2])
            };
            
            // Add to total field
            electricField[i][0] += E1[0];
            electricField[i][1] += E1[1];
            electricField[i][2] += E1[2];
            
            electricField[j][0] += E2[0];
            electricField[j][1] += E2[1];
            electricField[j][2] += E2[2];
        }
    }
}
*/  // End of DEPRECATED calculateInducedFieldS3S5

double DrudeSCFOM::updateDrudePositions(model::MCState& state,
                                            const std::vector<DrudeParticle>& particles,
                                            const std::vector<Vec3>& electricField,
                                            double damping,
                                            double maxStep,
                                            double hardWall) const {
    double maxDisplacement = 0.0;
    
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& particle = particles[i];
        auto& drude = state.atoms[particle.drudeIndex];
        const auto& parent = state.atoms[particle.parentIndex];
        
        // Freeze particles with zero polarizability
        if (particle.polarizability < 1e-14) {
            drude.x = parent.x;
            drude.y = parent.y;
            drude.z = parent.z;
            continue;
        }
        
        // Target displacement from force balance: F_electric = F_spring
        // The relationship μ = α * E_physical gives us the dipole moment
        // Since μ = |q_D| * d, we have d = α * E_physical / |q_D|
        // Our electric field includes ONE_4PI_EPS0, so E_physical = E_code / ONE_4PI_EPS0
        // Therefore: d = α * E_code / (|q_D| * ONE_4PI_EPS0)
        // But we also have k = q_D² * ONE_4PI_EPS0 / α
        // So: d = q_D * E_code / k gives the same result
        // This is correct! The issue must be elsewhere.
        double targetX = particle.charge * electricField[i][0] / particle.kSpring;
        double targetY = particle.charge * electricField[i][1] / particle.kSpring;
        double targetZ = particle.charge * electricField[i][2] / particle.kSpring;
        
        // Debug first particle with more detail
        static int debugIter = 0;
        if (i == 0 && debugIter < 5) {
            // Debug output with field magnitude information
            [[maybe_unused]] double fieldMag = std::sqrt(electricField[i][0]*electricField[i][0] + 
                                       electricField[i][1]*electricField[i][1] + 
                                       electricField[i][2]*electricField[i][2]);
            [[maybe_unused]] double targetMag = std::sqrt(targetX*targetX + targetY*targetY + targetZ*targetZ);
            
            // Show x-component details
            std::cout << "  UPDATE[0] iter " << debugIter 
                      << ": E_x=" << electricField[i][0]
                      << " target_dx=" << targetX
                      << " (q_D=" << particle.charge
                      << " k=" << particle.kSpring
                      << " parent.charge=" << parent.charge << ")"
                      << std::endl;
            debugIter++;
        }
        
        // Current displacement
        double oldX = drude.x - parent.x;
        double oldY = drude.y - parent.y;
        double oldZ = drude.z - parent.z;
        
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(oldX, oldY, oldZ, box);
        
        // Mixed update with damping
        double newX = (1.0 - damping) * oldX + damping * targetX;
        double newY = (1.0 - damping) * oldY + damping * targetY;
        double newZ = (1.0 - damping) * oldZ + damping * targetZ;
        
        // Apply step size limit
        double stepX = newX - oldX;
        double stepY = newY - oldY;
        double stepZ = newZ - oldZ;
        double step2 = stepX*stepX + stepY*stepY + stepZ*stepZ;
        
        if (step2 > maxStep * maxStep) {
            double scale = maxStep / std::sqrt(step2);
            newX = oldX + stepX * scale;
            newY = oldY + stepY * scale;
            newZ = oldZ + stepZ * scale;
        }
        
        // Apply hard wall constraint
        double r2 = newX*newX + newY*newY + newZ*newZ;
        if (hardWall < 1e8 && r2 > hardWall * hardWall) {
            double scale = hardWall / std::sqrt(r2);
            newX *= scale;
            newY *= scale;
            newZ *= scale;
        }
        
        // Update position
        drude.x = parent.x + newX;
        drude.y = parent.y + newY;
        drude.z = parent.z + newZ;
        
        // Track maximum displacement
        double change2 = (newX - oldX) * (newX - oldX) +
                        (newY - oldY) * (newY - oldY) +
                        (newZ - oldZ) * (newZ - oldZ);
        maxDisplacement = std::max(maxDisplacement, std::sqrt(change2));
    }
    
    return maxDisplacement;
}

double DrudeSCFOM::calculateSpringEnergy(const model::MCState& state,
                                             const std::vector<DrudeParticle>& particles) const {
    double energy = 0.0;
    
    for (const auto& p : particles) {
        const auto& drude = state.atoms[p.drudeIndex];
        const auto& parent = state.atoms[p.parentIndex];
        
        double dx = drude.x - parent.x;
        double dy = drude.y - parent.y;
        double dz = drude.z - parent.z;
        
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
        
        double r2 = dx*dx + dy*dy + dz*dz;
        energy += 0.5 * p.kSpring * r2;
    }
    
    return energy;
}

bool DrudeSCFOM::inSameMolecule(int atom1, int atom2, const model::MCState& state) const {
    // Check if atoms belong to same residue (molecule)
    // First check if we have residues at all
    if (state.residues.empty() || state.activeResidueCount == 0) {
        // No residue info means no intramolecular exclusions
        return false;
    }
    
    // Iterate through active residues
    for (int i = 0; i < state.activeResidueCount; ++i) {
        if (i >= static_cast<int>(state.residues.size())) {
            break;
        }
        
        const auto& res = state.residues[i];
        
        // Skip inactive residues (check both active flag and valid atomStart)
        if (!res.active && res.atomStart < 0) {
            continue;
        }
        
        int start = res.atomStart;
        int end = start + res.atomCount;
        
        // Check if both atoms are in this residue
        bool atom1InRes = (atom1 >= start && atom1 < end);
        bool atom2InRes = (atom2 >= start && atom2 < end);
        
        if (atom1InRes && atom2InRes) {
            return true;
        }
    }
    
    return false;
}

bool DrudeSCFOM::applyDIIS(DIISData& diis,
                           const std::vector<DrudeParticle>& particles,
                           const std::vector<Vec3>& electricField,
                           model::MCState& state) const {
    // DIIS (Direct Inversion in the Iterative Subspace) acceleration
    // Simplified implementation without external linear algebra library
    
    const int N = particles.size();
    if (N == 0) return false;
    
    // Collect current positions and compute residuals
    std::vector<double> currentPos(3 * N);
    std::vector<double> residual(3 * N);
    
    std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
    
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        const auto& drude = state.atoms[p.drudeIndex];
        const auto& parent = state.atoms[p.parentIndex];
        
        // Current displacement
        double dx = drude.x - parent.x;
        double dy = drude.y - parent.y;
        double dz = drude.z - parent.z;
        applyPBC(dx, dy, dz, box);
        
        currentPos[3*i] = dx;
        currentPos[3*i+1] = dy;
        currentPos[3*i+2] = dz;
        
        // Residual: r = q*E - k*d
        residual[3*i] = p.charge * electricField[i][0] - p.kSpring * dx;
        residual[3*i+1] = p.charge * electricField[i][1] - p.kSpring * dy;
        residual[3*i+2] = p.charge * electricField[i][2] - p.kSpring * dz;
    }
    
    // Add to history
    if (diis.historySize < diis.maxHistory) {
        diis.positions.push_back(currentPos);
        diis.residuals.push_back(residual);
        diis.historySize++;
    } else {
        // Replace oldest
        diis.positions.erase(diis.positions.begin());
        diis.residuals.erase(diis.residuals.begin());
        diis.positions.push_back(currentPos);
        diis.residuals.push_back(residual);
    }
    
    // Need at least 2 iterations for DIIS
    if (diis.historySize < 2 || !diis.enabled) {
        return false;
    }
    
    // For now, use a simplified Anderson mixing instead of full DIIS
    // This avoids the need for matrix inversion
    // Anderson mixing: x_{n+1} = (1-beta)*x_n + beta*G(x_n)
    // where G(x) is the fixed-point iteration
    
    const double mixingParam = 0.5;  // Anderson mixing parameter
    
    // Compute mixed position from last two iterations
    if (diis.historySize >= 2) {
        size_t last = diis.historySize - 1;
        size_t prev = diis.historySize - 2;
        
        // Compute residual norm for adaptive mixing
        double resNorm = 0.0;
        for (size_t j = 0; j < residual.size(); ++j) {
            resNorm += residual[j] * residual[j];
        }
        resNorm = std::sqrt(resNorm);
        
        // Adaptive mixing: use smaller mixing for larger residuals
        double adaptiveMixing = mixingParam;
        if (resNorm > 100.0) {
            adaptiveMixing *= 0.5;
        }
        
        // Anderson mixing update
        std::vector<double> newPos(3 * N);
        for (size_t j = 0; j < newPos.size(); ++j) {
            // Simple mixing between current and previous
            newPos[j] = (1.0 - adaptiveMixing) * diis.positions[last][j] + 
                        adaptiveMixing * diis.positions[prev][j];
        }
        
        // Update Drude positions
        for (size_t i = 0; i < particles.size(); ++i) {
            const auto& p = particles[i];
            auto& drude = state.atoms[p.drudeIndex];
            const auto& parent = state.atoms[p.parentIndex];
            
            drude.x = parent.x + newPos[3*i];
            drude.y = parent.y + newPos[3*i+1];
            drude.z = parent.z + newPos[3*i+2];
        }
        
        return true;
    }
    
    return false;
}

} // namespace exp
} // namespace cpu
} // namespace platform
} // namespace pygcmc