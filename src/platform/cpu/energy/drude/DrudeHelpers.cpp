// DrudeHelpers.cpp - Helper functions for Drude calculations

#include "DrudeForce.hpp"
#include "../common/EnergyUtils.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

void DrudeForce::updateDrudePositions(model::MCState& state, 
                                      const std::vector<Vec3>& displacements) {
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& particle = particles[i];
        state.atoms[particle.drudeIndex].x = state.atoms[particle.parentIndex].x + displacements[i].x;
        state.atoms[particle.drudeIndex].y = state.atoms[particle.parentIndex].y + displacements[i].y;
        state.atoms[particle.drudeIndex].z = state.atoms[particle.parentIndex].z + displacements[i].z;
    }
}

void DrudeForce::calculateElectricFieldResponse(model::MCState& state, 
                                               std::vector<Vec3>& response) {
    response.resize(particles.size());
    
    // Calculate electric field at each Drude
    std::vector<Vec3> electricField(particles.size(), Vec3(0.0, 0.0, 0.0));
    calculateElectricFieldAtDrudes(state, electricField, true);  // Include all charges
    
    // Convert field to displacement
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& particle = particles[i];
        double k = 4184.0 * 1000.0;  // Force constant
        
        response[i] = electricField[i] * (particle.charge * 138.935 / k);
    }
}

void DrudeForce::calculateElectricFieldAtDrudes(model::MCState& state,
                                               std::vector<Vec3>& electricField,
                                               bool includeDrudes) {
    // Calculate electric field at each Drude particle position
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& particle = particles[i];
        const auto& drude = state.atoms[particle.drudeIndex];
        
        Vec3 field(0.0, 0.0, 0.0);
        
        // Find which residue this Drude belongs to
        int drudeResIdx = -1;
        for (int resIdx = 0; resIdx < state.activeResidueCount; resIdx++) {
            const auto& res = state.residues[resIdx];
            if (particle.parentIndex >= res.atomStart && 
                particle.parentIndex < res.atomStart + res.atomCount) {
                drudeResIdx = resIdx;
                break;
            }
        }
        
        // Sum contributions from all charges
        for (int j = 0; j < state.activeAtomCount; j++) {
            const auto& atom = state.atoms[j];
            
            // Skip self
            if (j == particle.drudeIndex) continue;
            
            // Skip Drude particles if not included
            if (!includeDrudes) {
                bool isDrude = false;
                for (const auto& p : particles) {
                    if (j == p.drudeIndex) {
                        isDrude = true;
                        break;
                    }
                }
                if (isDrude) continue;
            }
            
            // Check intramolecular exclusion
            if (drudeResIdx >= 0) {
                const auto& res = state.residues[drudeResIdx];
                if (j >= res.atomStart && j < res.atomStart + res.atomCount) {
                    continue;  // Same residue - exclude
                }
            }
            
            // Calculate field contribution
            double dx = drude.x - atom.x;
            double dy = drude.y - atom.y;
            double dz = drude.z - atom.z;
            
            // Apply PBC
            applyPBC(dx, dy, dz, state.info.box[0], state.info.box[1], state.info.box[2]);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 < 0.01*0.01) continue;  // Too close
            
            double r = std::sqrt(r2);
            double factor = atom.charge / (r2 * r);
            
            field.x += factor * dx;
            field.y += factor * dy;
            field.z += factor * dz;
        }
        
        electricField[i] = field;
    }
}

double DrudeForce::calculateEnergyDirect(const model::MCState& state) {
    // Calculate total energy with current Drude positions
    // This is a simplified version - in production, use the full energy calculation
    
    double energy = 0.0;
    
    // Harmonic restraint energy
    for (const auto& particle : particles) {
        const auto& drude = state.atoms[particle.drudeIndex];
        const auto& parent = state.atoms[particle.parentIndex];
        
        double dx = drude.x - parent.x;
        double dy = drude.y - parent.y;
        double dz = drude.z - parent.z;
        
        double k = 4184.0 * 1000.0;  // kJ/mol/nm^2
        energy += 0.5 * k * (dx*dx + dy*dy + dz*dz);
    }
    
    // Electrostatic energy (simplified - doesn't include all terms)
    for (int i = 0; i < state.activeAtomCount; i++) {
        for (int j = i+1; j < state.activeAtomCount; j++) {
            const auto& atom1 = state.atoms[i];
            const auto& atom2 = state.atoms[j];
            
            double dx = atom1.x - atom2.x;
            double dy = atom1.y - atom2.y;
            double dz = atom1.z - atom2.z;
            
            applyPBC(dx, dy, dz, state.info.box[0], state.info.box[1], state.info.box[2]);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 > state.info.cutoff * state.info.cutoff) continue;
            
            double r = std::sqrt(r2);
            energy += 138.935 * atom1.charge * atom2.charge / r;
        }
    }
    
    return energy;
}

void DrudeForce::applyPBC(double& dx, double& dy, double& dz,
                          double boxX, double boxY, double boxZ) {
    // Apply periodic boundary conditions
    if (dx > boxX * 0.5) dx -= boxX;
    else if (dx < -boxX * 0.5) dx += boxX;
    
    if (dy > boxY * 0.5) dy -= boxY;
    else if (dy < -boxY * 0.5) dy += boxY;
    
    if (dz > boxZ * 0.5) dz -= boxZ;
    else if (dz < -boxZ * 0.5) dz += boxZ;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc