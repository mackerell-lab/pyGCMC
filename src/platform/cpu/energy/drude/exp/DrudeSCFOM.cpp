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
        m_lastIterationCount = 0;
        return true;
    }
    
    // Reset iteration counter
    m_lastIterationCount = 0;
    
    // Initialize electric field vectors
    std::vector<Vec3> electricField(particles.size(), {0.0, 0.0, 0.0});
    
    // DIIS acceleration data
    DIISData diis;
    diis.maxHistory = 8;  // Default max history
    
    // History for residuals and displacements
    std::vector<double> lastDisp(3*particles.size(), 0.0);
    double prevMaxRes = 1e300;
    
    // Estimate spectral radius for adaptive damping
    double rho = estimateSpectralRadius(state, particles, pairs);
    
    // Adaptive damping parameter
    double dampingFactor = params.dampingFactor;
    if (params.enableAdaptiveDamping) {
        dampingFactor = computeAdaptiveDamping(rho);
        if (params.logLevel >= 1) {
            std::cout << "Spectral radius ρ = " << rho 
                     << ", using adaptive damping = " << dampingFactor << std::endl;
        }
    }
    
    // Get initial positions
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        const auto& drude = state.atoms[p.drudeIndex];
        const auto& parent = state.atoms[p.parentIndex];
        double dx = drude.x - parent.x;
        double dy = drude.y - parent.y;
        double dz = drude.z - parent.z;
        std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
        applyPBC(dx, dy, dz, box);
        lastDisp[3*i] = dx;
        lastDisp[3*i+1] = dy;
        lastDisp[3*i+2] = dz;
    }
    
    // SCF iteration
    for (int iter = 0; iter < params.maxIterations; ++iter) {
        m_lastIterationCount = iter + 1;  // Update iteration count
        
        // 1. Calculate total electric field (unified field including all contributions)
        std::fill(electricField.begin(), electricField.end(), Vec3{0.0, 0.0, 0.0});
        calculateElectricField(state, particles, pairs, electricField, params);
        
        // 2. Calculate residuals and displacement changes
        double maxRes = 0.0, maxDelta = 0.0;
        std::vector<double> currentDisp(3*particles.size());
        std::vector<double> currentRes(3*particles.size());
        
        for (size_t i = 0; i < particles.size(); ++i) {
            const auto& p = particles[i];
            if (p.polarizability < 1e-14) {
                // Frozen particle
                currentDisp[3*i] = 0.0;
                currentDisp[3*i+1] = 0.0;
                currentDisp[3*i+2] = 0.0;
                currentRes[3*i] = 0.0;
                currentRes[3*i+1] = 0.0;
                currentRes[3*i+2] = 0.0;
                continue;
            }
            
            const auto& drude = state.atoms[p.drudeIndex];
            const auto& parent = state.atoms[p.parentIndex];
            
            double dx = drude.x - parent.x;
            double dy = drude.y - parent.y;
            double dz = drude.z - parent.z;
            
            std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
            applyPBC(dx, dy, dz, box);
            
            // Residual: r = qE - kd
            double rx = p.charge * electricField[i][0] - p.kSpring * dx;
            double ry = p.charge * electricField[i][1] - p.kSpring * dy;
            double rz = p.charge * electricField[i][2] - p.kSpring * dz;
            double resNorm = std::sqrt(rx*rx + ry*ry + rz*rz);
            maxRes = std::max(maxRes, resNorm);
            
            // Displacement change
            double ddx = dx - lastDisp[3*i];
            double ddy = dy - lastDisp[3*i+1];
            double ddz = dz - lastDisp[3*i+2];
            double deltaNorm = std::sqrt(ddx*ddx + ddy*ddy + ddz*ddz);
            maxDelta = std::max(maxDelta, deltaNorm);
            
            currentDisp[3*i] = dx;
            currentDisp[3*i+1] = dy;
            currentDisp[3*i+2] = dz;
            currentRes[3*i] = rx;
            currentRes[3*i+1] = ry;
            currentRes[3*i+2] = rz;
        }
        
        // 3. Check convergence (residual-based)
        if (maxRes < params.tolerance) {
            if (params.logLevel > 0) {
                std::cout << "SCF converged at iter " << iter 
                         << ", maxRes = " << maxRes << std::endl;
            }
            return true;
        }
        
        // Alternative convergence: displacement change
        if (maxDelta < params.displacementTolerance && iter > 5) {
            if (params.logLevel > 0) {
                std::cout << "SCF converged (displacement) at iter " << iter 
                         << ", maxDelta = " << maxDelta << std::endl;
            }
            return true;
        }
        
        // 4. DIIS acceleration (if conditions are met)
        bool diisApplied = false;
        if (iter >= params.diisStartIter && maxRes < prevMaxRes) {
            diis.positions.push_back(currentDisp);
            diis.residuals.push_back(currentRes);
            if (diis.positions.size() > (size_t)params.diisMaxHistory) {
                diis.positions.erase(diis.positions.begin());
                diis.residuals.erase(diis.residuals.begin());
            }
            diisApplied = applyDIIS(diis, particles, state);
        }
        
        // 5. If DIIS not applied, use adaptive damping with backtracking
        if (!diisApplied) {
            double lambda = dampingFactor;
            
            // Save current positions
            std::vector<Vec3> savedPos(particles.size());
            for (size_t i = 0; i < particles.size(); ++i) {
                const auto& p = particles[i];
                const auto& drude = state.atoms[p.drudeIndex];
                savedPos[i] = {drude.x, drude.y, drude.z};
            }
            
            // Backtracking line search
            for (int ls = 0; ls < 6; ++ls) {
                // Update positions
                for (size_t i = 0; i < particles.size(); ++i) {
                    const auto& p = particles[i];
                    if (p.polarizability < 1e-14) continue;
                    
                    auto& drude = state.atoms[p.drudeIndex];
                    const auto& parent = state.atoms[p.parentIndex];
                    
                    // Target displacement
                    // k already contains ONE_4PI_EPS0, so we don't need to multiply it again
                    // d = q_D * E / k
                    double targetX = p.charge * electricField[i][0] / p.kSpring;
                    double targetY = p.charge * electricField[i][1] / p.kSpring;
                    double targetZ = p.charge * electricField[i][2] / p.kSpring;
                    
                    // Mixed update
                    double newX = (1.0 - lambda) * currentDisp[3*i] + lambda * targetX;
                    double newY = (1.0 - lambda) * currentDisp[3*i+1] + lambda * targetY;
                    double newZ = (1.0 - lambda) * currentDisp[3*i+2] + lambda * targetZ;
                    
                    // Limit step size
                    double step2 = (newX-currentDisp[3*i])*(newX-currentDisp[3*i]) +
                                  (newY-currentDisp[3*i+1])*(newY-currentDisp[3*i+1]) +
                                  (newZ-currentDisp[3*i+2])*(newZ-currentDisp[3*i+2]);
                    if (step2 > params.maxStep * params.maxStep) {
                        double scale = params.maxStep / std::sqrt(step2);
                        newX = currentDisp[3*i] + scale * (newX - currentDisp[3*i]);
                        newY = currentDisp[3*i+1] + scale * (newY - currentDisp[3*i+1]);
                        newZ = currentDisp[3*i+2] + scale * (newZ - currentDisp[3*i+2]);
                    }
                    
                    // Apply hard wall if enabled
                    if (params.enableHardWall) {
                        double r2 = newX*newX + newY*newY + newZ*newZ;
                        if (r2 > params.maxDrudeDistance * params.maxDrudeDistance) {
                            double scale = params.maxDrudeDistance / std::sqrt(r2);
                            newX *= scale;
                            newY *= scale;
                            newZ *= scale;
                        }
                    }
                    
                    drude.x = parent.x + newX;
                    drude.y = parent.y + newY;
                    drude.z = parent.z + newZ;
                }
                
                // Recalculate field and residual
                std::fill(electricField.begin(), electricField.end(), Vec3{0.0, 0.0, 0.0});
                calculateElectricField(state, particles, pairs, electricField, params);
                
                double newMaxRes = 0.0;
                for (size_t i = 0; i < particles.size(); ++i) {
                    const auto& p = particles[i];
                    if (p.polarizability < 1e-14) continue;
                    
                    const auto& drude = state.atoms[p.drudeIndex];
                    const auto& parent = state.atoms[p.parentIndex];
                    
                    double dx = drude.x - parent.x;
                    double dy = drude.y - parent.y;
                    double dz = drude.z - parent.z;
                    
                    std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
                    applyPBC(dx, dy, dz, box);
                    
                    double rx = p.charge * electricField[i][0] - p.kSpring * dx;
                    double ry = p.charge * electricField[i][1] - p.kSpring * dy;
                    double rz = p.charge * electricField[i][2] - p.kSpring * dz;
                    double resNorm = std::sqrt(rx*rx + ry*ry + rz*rz);
                    newMaxRes = std::max(newMaxRes, resNorm);
                }
                
                if (newMaxRes < 0.98 * maxRes || ls == 5) {
                    maxRes = newMaxRes;
                    if (ls == 0 && newMaxRes < 0.5 * maxRes) {
                        dampingFactor = std::min(0.8, lambda * 1.1);
                    }
                    break;
                } else {
                    lambda *= 0.5;
                    // Restore positions
                    for (size_t i = 0; i < particles.size(); ++i) {
                        const auto& p = particles[i];
                        auto& drude = state.atoms[p.drudeIndex];
                        drude.x = savedPos[i][0];
                        drude.y = savedPos[i][1];
                        drude.z = savedPos[i][2];
                    }
                }
            }
        }
        
        // 6. Update history
        lastDisp = currentDisp;
        prevMaxRes = maxRes;
        
        // 7. Logging
        if (params.logLevel >= 2) {
            std::cout << "Iter " << iter << ": maxRes = " << maxRes 
                     << ", maxDelta = " << maxDelta << std::endl;
        }
    }
    
    if (params.logLevel > 0) {
        std::cout << "SCF did not converge after " << params.maxIterations 
                 << " iterations, final maxRes = " << prevMaxRes << std::endl;
    }
    return false;  // Did not converge
}

void DrudeSCFOM::calculateElectricField(const model::MCState& state,
                                            const std::vector<DrudeParticle>& particles,
                                            const std::vector<ScreenedPair>& pairs,
                                            std::vector<Vec3>& electricField,
                                            const DrudeSCFParams& params) const {
    // Calculate total electric field at each Drude particle
    // Field = External field (from non-Drude charges) + Induced field (from other dipoles)
    
    // Clear field
    for (auto& field : electricField) {
        field[0] = field[1] = field[2] = 0.0;
    }
    
    // Build maps for quick lookup
    std::map<std::pair<int,int>, double> screenedPairs;
    std::vector<bool> isDrude(state.activeAtomCount, false);
    std::vector<int> atomToDipole(state.activeAtomCount, -1);
    
    // Mark Drude atoms and build mapping
    for (size_t i = 0; i < particles.size(); ++i) {
        isDrude[particles[i].drudeIndex] = true;
        atomToDipole[particles[i].parentIndex] = i;
        atomToDipole[particles[i].drudeIndex] = i;
    }
    
    // Build screened pairs map
    for (const auto& pair : pairs) {
        screenedPairs[{pair.dipole1, pair.dipole2}] = pair.thole;
        screenedPairs[{pair.dipole2, pair.dipole1}] = pair.thole;
    }
    
    std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
    
    // Calculate field at each PARENT position (not Drude position!)
    // This is critical for correct SCF convergence
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& pi = particles[i];
        const auto& fieldPoint = state.atoms[pi.parentIndex];  // Field at PARENT position
        const auto& parentI = state.atoms[pi.parentIndex];
        
        Vec3 field = {0.0, 0.0, 0.0};
        
        // 1. External field from non-Drude charges
        double extFieldX = 0, extFieldY = 0, extFieldZ = 0;
        for (int j = 0; j < state.activeAtomCount; ++j) {
            // Skip self and Drude charges
            if (j == pi.parentIndex || j == pi.drudeIndex || isDrude[j]) continue;
            
            const auto& atomJ = state.atoms[j];
            if (std::abs(atomJ.charge) < 1e-14) continue;
            
            if (i == 0 && params.logLevel >= 3) {
                std::cout << "  External source j=" << j 
                         << " charge=" << atomJ.charge 
                         << " isDrude=" << isDrude[j] << std::endl;
            }
            
            // Determine screening factor
            double S1 = 1.0;  // Default: no screening
            
            // Check if j belongs to a dipole that has screening with dipole i
            int dipoleK = atomToDipole[j];
            if (dipoleK >= 0 && dipoleK < static_cast<int>(particles.size())) {
                auto it = screenedPairs.find({static_cast<int>(i), dipoleK});
                if (it != screenedPairs.end()) {
                    // Compute S1 using parent-parent distance
                    const auto& parentK = state.atoms[particles[dipoleK].parentIndex];
                    
                    double dx_pp = parentK.x - parentI.x;
                    double dy_pp = parentK.y - parentI.y;
                    double dz_pp = parentK.z - parentI.z;
                    applyPBC(dx_pp, dy_pp, dz_pp, box);
                    double r_pp = std::sqrt(dx_pp*dx_pp + dy_pp*dy_pp + dz_pp*dz_pp);
                    
                    if (r_pp > 1e-12) {
                        S1 = tholeS1(r_pp, pi.polarizability, 
                                    particles[dipoleK].polarizability, it->second);
                    }
                }
            }
            
            // Calculate field contribution (field_point - source)
            double dx = fieldPoint.x - atomJ.x;
            double dy = fieldPoint.y - atomJ.y;
            double dz = fieldPoint.z - atomJ.z;
            applyPBC(dx, dy, dz, box);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1e-12) continue;
            
            double r = std::sqrt(r2);
            double r3 = r2 * r;
            
            // E = k * q * S(u) / r^2 * r_hat
            double factor = DrudeConstants::ONE_4PI_EPS0 * atomJ.charge * S1 / r3;
            
            extFieldX += factor * dx;
            extFieldY += factor * dy;
            extFieldZ += factor * dz;
        }
        
        field[0] += extFieldX;
        field[1] += extFieldY;
        field[2] += extFieldZ;
        
        if (i == 0 && params.logLevel >= 3) {
            std::cout << "  External field: (" << extFieldX << ", " 
                     << extFieldY << ", " << extFieldZ << ")" << std::endl;
        }
        
        // 2. Induced field from OTHER dipoles (using S1 point charge model)
        for (size_t j = 0; j < particles.size(); ++j) {
            if (i == j) continue;  // Skip self
            
            const auto& pj = particles[j];
            const auto& drudeJ = state.atoms[pj.drudeIndex];
            const auto& parentJ = state.atoms[pj.parentIndex];
            
            // Get screening parameter if this is a screened pair
            double thole = 0.0;
            auto it = screenedPairs.find({static_cast<int>(i), static_cast<int>(j)});
            if (it != screenedPairs.end()) {
                thole = it->second;
            }
            
            // Calculate parent-parent distance for screening
            double dx_pp = parentJ.x - parentI.x;
            double dy_pp = parentJ.y - parentI.y;
            double dz_pp = parentJ.z - parentI.z;
            applyPBC(dx_pp, dy_pp, dz_pp, box);
            double r_pp = std::sqrt(dx_pp*dx_pp + dy_pp*dy_pp + dz_pp*dz_pp);
            
            if (r_pp < 1e-12) continue;
            
            // Calculate S1 screening for all 4 interactions
            double S1 = 1.0;
            if (thole > 1e-10) {
                S1 = tholeS1(r_pp, pi.polarizability, pj.polarizability, thole);
            }
            
            // Field from Drude charge of dipole j
            double dx_dj = fieldPoint.x - drudeJ.x;
            double dy_dj = fieldPoint.y - drudeJ.y;
            double dz_dj = fieldPoint.z - drudeJ.z;
            applyPBC(dx_dj, dy_dj, dz_dj, box);
            
            double r2_dj = dx_dj*dx_dj + dy_dj*dy_dj + dz_dj*dz_dj;
            if (r2_dj > 1e-12) {
                double r_dj = std::sqrt(r2_dj);
                double r3_dj = r2_dj * r_dj;
                double factor_dj = DrudeConstants::ONE_4PI_EPS0 * pj.charge * S1 / r3_dj;
                field[0] += factor_dj * dx_dj;
                field[1] += factor_dj * dy_dj;
                field[2] += factor_dj * dz_dj;
            }
            
            // Parent charges are already included in external field, don't double count!
            // Only include Drude charges in induced field
        }
        
        electricField[i] = field;
        
        // Debug output
        if (i == 0 && params.logLevel >= 2) {
            double mag = std::sqrt(field[0]*field[0] + field[1]*field[1] + field[2]*field[2]);
            std::cout << "TOTAL_FIELD[0]: E_x=" << field[0] 
                     << " E_y=" << field[1] 
                     << " E_z=" << field[2]
                     << " |E|=" << mag << std::endl;
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
        // The physical relationship: μ = α * E gives the dipole moment
        // Since μ = |q_D| * d, we have d = α * E / |q_D|
        // In MD units with k = q_D² * ONE_4PI_EPS0 / α:
        // k already contains ONE_4PI_EPS0, so d = q_D * E / k
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
                           model::MCState& state) const {
    // Real DIIS (Direct Inversion in the Iterative Subspace) acceleration
    // Based on Pulay's method using residual dot products
    
    const int N = particles.size();
    if (N == 0 || diis.positions.size() < 3) return false;
    
    int m = diis.positions.size();
    
    // Construct B matrix: B_ij = <r_i, r_j>
    std::vector<double> B((m+1)*(m+1), 0.0);
    std::vector<double> rhs(m+1, 0.0);
    
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            double dotProduct = 0.0;
            for (size_t k = 0; k < diis.residuals[i].size(); ++k) {
                dotProduct += diis.residuals[i][k] * diis.residuals[j][k];
            }
            B[i*(m+1)+j] = dotProduct;
        }
        B[i*(m+1)+m] = -1.0;
        B[m*(m+1)+i] = -1.0;
    }
    rhs[m] = -1.0;
    
    // Solve linear system using Gaussian elimination
    std::vector<double> c(m+1);
    
    // Simple Gaussian elimination (without pivoting for simplicity)
    for (int i = 0; i < m+1; ++i) {
        // Find pivot
        double pivot = B[i*(m+1)+i];
        if (std::abs(pivot) < 1e-10) {
            // Matrix is singular, fall back to no DIIS
            return false;
        }
        
        // Scale row
        for (int j = i; j < m+1; ++j) {
            B[i*(m+1)+j] /= pivot;
        }
        rhs[i] /= pivot;
        
        // Eliminate column
        for (int k = i+1; k < m+1; ++k) {
            double factor = B[k*(m+1)+i];
            for (int j = i; j < m+1; ++j) {
                B[k*(m+1)+j] -= factor * B[i*(m+1)+j];
            }
            rhs[k] -= factor * rhs[i];
        }
    }
    
    // Back substitution
    for (int i = m; i >= 0; --i) {
        c[i] = rhs[i];
        for (int j = i+1; j < m+1; ++j) {
            c[i] -= B[i*(m+1)+j] * c[j];
        }
    }
    
    // Mix positions: x_new = Σ c_i * x_i
    std::vector<double> newPos(3*N, 0.0);
    for (int i = 0; i < m; ++i) {
        for (size_t k = 0; k < newPos.size(); ++k) {
            newPos[k] += c[i] * diis.positions[i][k];
        }
    }
    
    // Update Drude positions
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        if (p.polarizability < 1e-14) continue;
        
        auto& drude = state.atoms[p.drudeIndex];
        const auto& parent = state.atoms[p.parentIndex];
        
        // Apply max step constraint
        double dx = newPos[3*i];
        double dy = newPos[3*i+1];
        double dz = newPos[3*i+2];
        
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 > 0.04) {  // Max 0.2 nm
            double scale = 0.2 / std::sqrt(r2);
            dx *= scale;
            dy *= scale;
            dz *= scale;
        }
        
        drude.x = parent.x + dx;
        drude.y = parent.y + dy;
        drude.z = parent.z + dz;
    }
    
    return true;
}

double DrudeSCFOM::estimateSpectralRadius(const model::MCState& state,
                                          const std::vector<DrudeParticle>& particles,
                                          const std::vector<ScreenedPair>& pairs) const {
    // Estimate spectral radius ρ(αT) for the system
    // For simplicity, use pairwise approximation: ρ ≈ max(α_i * T_ij)
    
    if (particles.size() < 2) {
        return 0.0;  // Single dipole has no mutual polarization
    }
    
    double maxRho = 0.0;
    std::array<double, 3> box = {state.info.box[0], state.info.box[1], state.info.box[2]};
    
    // Build screened pairs map
    std::map<std::pair<int,int>, double> screenedPairs;
    for (const auto& pair : pairs) {
        screenedPairs[{pair.dipole1, pair.dipole2}] = pair.thole;
        screenedPairs[{pair.dipole2, pair.dipole1}] = pair.thole;
    }
    
    // Check all dipole pairs
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& pi = particles[i];
        if (pi.polarizability < 1e-14) continue;
        
        const auto& parent1 = state.atoms[pi.parentIndex];
        
        for (size_t j = 0; j < particles.size(); ++j) {
            if (i == j) continue;
            
            const auto& pj = particles[j];
            if (pj.polarizability < 1e-14) continue;
            
            const auto& parent2 = state.atoms[pj.parentIndex];
            
            // Calculate distance between parents
            double dx = parent2.x - parent1.x;
            double dy = parent2.y - parent1.y;
            double dz = parent2.z - parent1.z;
            applyPBC(dx, dy, dz, box);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 1e-12) continue;
            
            double r = std::sqrt(r2);
            
            // Get screening parameter
            double a_pair = 0.0;
            auto it = screenedPairs.find({static_cast<int>(i), static_cast<int>(j)});
            if (it != screenedPairs.end()) {
                a_pair = it->second;
            }
            
            // Calculate screening
            double s = tholeS1(r, pi.polarizability, pj.polarizability, a_pair);
            
            // Dipole-dipole interaction strength
            double T = DrudeConstants::ONE_4PI_EPS0 * s / (r*r*r);
            
            // Spectral radius contribution
            double rho_ij = std::sqrt(pi.polarizability * pj.polarizability) * T;
            maxRho = std::max(maxRho, rho_ij);
        }
    }
    
    return maxRho;
}

double DrudeSCFOM::computeAdaptiveDamping(double rho) const {
    // Compute adaptive damping factor based on spectral radius
    // Ensures stable convergence even near instability
    
    if (rho < 0.9) {
        // Stable region: standard damping
        return 0.5;
    } else if (rho < 1.0) {
        // Near instability: scale damping to maintain stability
        // damping = 0.9 / rho ensures effective ρ < 0.9
        return std::min(0.5, 0.9 / rho);
    } else if (rho < 2.0) {
        // Unstable but tractable: strong damping
        return std::min(0.3, 0.9 / rho);
    } else {
        // Highly unstable: very strong damping
        return std::min(0.1, 0.5 / rho);
    }
}

} // namespace exp
} // namespace cpu
} // namespace platform
} // namespace pygcmc