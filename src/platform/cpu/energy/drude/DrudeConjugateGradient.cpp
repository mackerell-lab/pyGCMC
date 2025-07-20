#include "DrudeConjugateGradient.hpp"
#include <iostream>
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

using Vec3 = DrudeForce::Vec3;

bool DrudeConjugateGradient::minimizeDrudePositions(
    model::MCState& state,
    const DrudeForce& drudeForce,
    double tolerance,
    int maxIterations) {
    
    const auto& particles = drudeForce.particles;
    int nDrudes = particles.size();
    
    if (nDrudes == 0) return true;
    
    // Initialize Drude positions at parent positions
    for (const auto& drude : particles) {
        state.atoms[drude.drudeIndex].x = state.atoms[drude.parentIndex].x;
        state.atoms[drude.drudeIndex].y = state.atoms[drude.parentIndex].y;
        state.atoms[drude.drudeIndex].z = state.atoms[drude.parentIndex].z;
    }
    
    // Initial solution vector (displacement from parent)
    std::vector<Vec3> x(nDrudes, Vec3(0, 0, 0));
    
    // Calculate right-hand side: b = F_external
    std::vector<Vec3> b = calculateRHS(state, drudeForce);
    
    // Initial residual: r = b - A*x = b (since x=0)
    std::vector<Vec3> r = b;
    
    // Initial search direction: p = r
    std::vector<Vec3> p = r;
    
    // Initial residual norm squared
    double rsold = dotProduct(r, r);
    
    // Check if already converged
    double initialNorm = std::sqrt(rsold / (3.0 * nDrudes));
    if (initialNorm < tolerance) {
        return true;
    }
    
    // Conjugate gradient iterations
    for (int iter = 0; iter < maxIterations; ++iter) {
        // Apply system matrix: Ap = A * p
        std::vector<Vec3> Ap = applySystemMatrix(state, drudeForce, p);
        
        // Calculate step size: alpha = r^T*r / p^T*A*p
        double pAp = dotProduct(p, Ap);
        
        // Check for breakdown
        if (std::abs(pAp) < 1e-10) {
            std::cerr << "CG breakdown: p^T*A*p near zero" << std::endl;
            break;
        }
        
        double alpha = rsold / pAp;
        
        // Update solution: x = x + alpha * p
        for (int i = 0; i < nDrudes; ++i) {
            x[i] += p[i] * alpha;
        }
        
        // Update residual: r = r - alpha * A*p
        for (int i = 0; i < nDrudes; ++i) {
            r[i] -= Ap[i] * alpha;
        }
        
        // New residual norm squared
        double rsnew = dotProduct(r, r);
        
        // Check convergence
        double rmsResidual = std::sqrt(rsnew / (3.0 * nDrudes));
        if (rmsResidual < tolerance) {
            // Update positions and verify convergence
            updatePositions(state, drudeForce, x);
            
            // Final force check
            double maxForce = calculateMaxForce(state, drudeForce);
            if (maxForce < tolerance) {
                return true;
            }
        }
        
        // Update search direction: p = r + beta * p
        double beta = rsnew / rsold;
        for (int i = 0; i < nDrudes; ++i) {
            p[i] = r[i] + p[i] * beta;
        }
        
        rsold = rsnew;
        
        // Apply hard wall constraint during iterations
        for (int i = 0; i < nDrudes; ++i) {
            const auto& drude = particles[i];
            double displacement = x[i].norm();
            
            if (displacement > drudeForce.scfParams.maxDrudeDistance) {
                double scale = drudeForce.scfParams.maxDrudeDistance / displacement;
                x[i] = x[i] * scale;
                
                // Reset CG after constraint (optional)
                r = calculateRHS(state, drudeForce);
                std::vector<Vec3> Ax = applySystemMatrix(state, drudeForce, x);
                for (int j = 0; j < nDrudes; ++j) {
                    r[j] = r[j] - Ax[j];
                }
                p = r;
                rsold = dotProduct(r, r);
            }
        }
    }
    
    // Update positions even if not fully converged
    updatePositions(state, drudeForce, x);
    
    return false;
}

std::vector<DrudeForce::Vec3> DrudeConjugateGradient::applySystemMatrix(
    const model::MCState& state,
    const DrudeForce& drudeForce,
    const std::vector<DrudeForce::Vec3>& x) {
    
    const auto& particles = drudeForce.particles;
    const auto& screenedPairs = drudeForce.screenedPairs;
    int nDrudes = particles.size();
    
    std::vector<Vec3> result(nDrudes, Vec3(0, 0, 0));
    
    // 1. Spring force contribution: K_spring * x
    for (int i = 0; i < nDrudes; ++i) {
        const auto& drude = particles[i];
        
        // Isotropic spring constant
        result[i] = x[i] * drude.kIsotropic;
        
        // TODO: Add anisotropic contributions if needed
    }
    
    // 2. Thole interaction contributions
    const double COULOMB_CONSTANT = 138.935456; // kJ*nm/mol/e^2
    
    for (const auto& pair : screenedPairs) {
        int i = pair.dipole1;
        int j = pair.dipole2;
        
        const auto& drude1 = particles[i];
        const auto& drude2 = particles[j];
        
        // Parent positions
        Vec3 pos1(state.atoms[drude1.parentIndex].x,
                  state.atoms[drude1.parentIndex].y,
                  state.atoms[drude1.parentIndex].z);
        Vec3 pos2(state.atoms[drude2.parentIndex].x,
                  state.atoms[drude2.parentIndex].y,
                  state.atoms[drude2.parentIndex].z);
        
        // Distance between parents
        Vec3 r12 = pos2 - pos1;
        
        // Apply periodic boundary conditions
        if (state.info.box[0] > 0) {
            r12.x -= state.info.box[0] * std::round(r12.x / state.info.box[0]);
            r12.y -= state.info.box[1] * std::round(r12.y / state.info.box[1]);
            r12.z -= state.info.box[2] * std::round(r12.z / state.info.box[2]);
        }
        
        double r = r12.norm();
        if (r < state.info.cutoff && r > 1e-6) {
            // Thole screening parameter
            double alpha1 = drude1.polarizability;
            double alpha2 = drude2.polarizability;
            double uscale = pair.thole / std::pow(alpha1 * alpha2, 1.0/6.0);
            double u = r * uscale;
            double u3 = u * u * u;
            
            // Thole damping function and derivatives
            double screening = 1.0 - std::exp(-u3);
            double dscreening = 3.0 * u * u * uscale * std::exp(-u3);
            
            // Dipole-dipole interaction tensor
            double r2 = r * r;
            double r3 = r2 * r;
            double r5 = r3 * r2;
            
            // T_ij = (3*r_i*r_j/r^5 - delta_ij/r^3) * screening
            // Apply T_ij to displacement j and add to force on i
            Vec3 rhat = r12 * (1.0 / r);
            
            // Interaction between induced dipoles
            double q1 = drude1.charge;
            double q2 = drude2.charge;
            
            // Matrix element contribution
            Vec3 Tij_xj = (rhat * (3.0 * rhat.dot(x[j]) / r5) - x[j] * (1.0 / r3)) * screening;
            result[i] -= Tij_xj * (COULOMB_CONSTANT * q1 * q2);
            
            // By symmetry
            Vec3 Tji_xi = (rhat * (3.0 * rhat.dot(x[i]) / r5) - x[i] * (1.0 / r3)) * screening;
            result[j] -= Tji_xi * (COULOMB_CONSTANT * q1 * q2);
        }
    }
    
    return result;
}

std::vector<DrudeForce::Vec3> DrudeConjugateGradient::calculateRHS(
    const model::MCState& state,
    const DrudeForce& drudeForce) {
    
    const auto& particles = drudeForce.particles;
    int nDrudes = particles.size();
    int nAtoms = state.activeAtomCount;
    
    std::vector<Vec3> rhs(nDrudes, Vec3(0, 0, 0));
    
    const double COULOMB_CONSTANT = 138.935456;
    
    // Calculate electric field at each Drude position from all other atoms
    // (excluding atoms in the same residue)
    for (int i = 0; i < nDrudes; ++i) {
        const auto& drude = particles[i];
        int drudeIdx = drude.drudeIndex;
        int parentIdx = drude.parentIndex;
        double drudeCharge = drude.charge;
        
        // Find which residue this Drude belongs to
        int drudeResidueIdx = -1;
        for (int resIdx = 0; resIdx < state.activeResidueCount; ++resIdx) {
            const auto& res = state.residues[resIdx];
            if (drudeIdx >= res.atomStart && drudeIdx < res.atomStart + res.atomCount) {
                drudeResidueIdx = resIdx;
                break;
            }
        }
        
        // Drude position (initially at parent)
        Vec3 drudePos(state.atoms[parentIdx].x,
                      state.atoms[parentIdx].y,
                      state.atoms[parentIdx].z);
        
        Vec3 electricField(0, 0, 0);
        
        // Sum electric field from all atoms NOT in the same residue
        for (int j = 0; j < nAtoms; ++j) {
            // Skip if atom is in the same residue
            bool sameResidue = false;
            if (drudeResidueIdx >= 0) {
                const auto& res = state.residues[drudeResidueIdx];
                if (j >= res.atomStart && j < res.atomStart + res.atomCount) {
                    sameResidue = true;
                }
            }
            if (sameResidue) continue;
            
            Vec3 atomPos(state.atoms[j].x,
                        state.atoms[j].y,
                        state.atoms[j].z);
            double atomCharge = state.atoms[j].charge;
            
            if (std::abs(atomCharge) < 1e-10) continue;
            
            Vec3 r = atomPos - drudePos;
            
            // Periodic boundary conditions
            if (state.info.box[0] > 0) {
                r.x -= state.info.box[0] * std::round(r.x / state.info.box[0]);
                r.y -= state.info.box[1] * std::round(r.y / state.info.box[1]);
                r.z -= state.info.box[2] * std::round(r.z / state.info.box[2]);
            }
            
            double dist = r.norm();
            if (dist < state.info.cutoff && dist > 1e-6) {
                // E = k * q / r^2 * r_hat
                electricField += r * (COULOMB_CONSTANT * atomCharge / (dist * dist * dist));
            }
        }
        
        // Force = charge * field
        rhs[i] = electricField * drudeCharge;
    }
    
    return rhs;
}

double DrudeConjugateGradient::dotProduct(
    const std::vector<DrudeForce::Vec3>& a,
    const std::vector<DrudeForce::Vec3>& b) {
    
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += a[i].dot(b[i]);
    }
    return sum;
}

void DrudeConjugateGradient::updatePositions(
    model::MCState& state,
    const DrudeForce& drudeForce,
    const std::vector<DrudeForce::Vec3>& solution) {
    
    const auto& particles = drudeForce.particles;
    
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& drude = particles[i];
        int drudeIdx = drude.drudeIndex;
        int parentIdx = drude.parentIndex;
        
        // Update Drude position: r_drude = r_parent + delta_r
        state.atoms[drudeIdx].x = state.atoms[parentIdx].x + solution[i].x;
        state.atoms[drudeIdx].y = state.atoms[parentIdx].y + solution[i].y;
        state.atoms[drudeIdx].z = state.atoms[parentIdx].z + solution[i].z;
    }
}

double DrudeConjugateGradient::calculateMaxForce(
    const model::MCState& state,
    const DrudeForce& drudeForce) {
    
    // Calculate forces on all atoms including Drudes
    std::vector<Vec3> forces(state.activeAtomCount, Vec3(0, 0, 0));
    
    // Use DrudeForce's own force calculation
    // Create a mutable copy of state for force calculation
    model::MCState stateCopy = state;
    const_cast<DrudeForce&>(drudeForce).calculateForces(stateCopy, forces);
    
    // Find maximum force on Drude particles
    double maxForce = 0.0;
    const auto& particles = drudeForce.particles;
    
    for (const auto& drude : particles) {
        double forceMag = forces[drude.drudeIndex].norm();
        maxForce = std::max(maxForce, forceMag);
    }
    
    return maxForce;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc