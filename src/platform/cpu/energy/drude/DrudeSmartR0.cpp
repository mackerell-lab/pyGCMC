// DrudeSmartR0.cpp - Smart r0 correction for improved OPT3 accuracy

#include "DrudeForce.hpp"
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * Calculate "smart" r0 that includes statistical correction for 
 * missing Drude-Drude interactions
 */
void DrudeForce::calculateSmartR0(model::MCState& state, std::vector<Vec3>& r0) {
    // First calculate standard r0 (static field from non-Drude atoms)
    calculateStandardR0(state, r0);
    
    // Now apply corrections based on local environment
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& particle = particles[i];
        const auto& drude = state.atoms[particle.drudeIndex];
        
        // 1. Calculate local density metric
        double localDensity = calculateLocalDensity(state, drude, i);
        
        // 2. Check for hydrogen bonding environment
        bool inHBondNetwork = isInHydrogenBondNetwork(state, particle);
        
        // 3. Calculate charge asymmetry in local environment
        double chargeAsymmetry = calculateLocalChargeAsymmetry(state, drude, i);
        
        // 4. Apply correction factors
        double correctionFactor = 1.0;
        
        // Density correction: more neighbors = stronger Drude-Drude interactions
        // Empirical formula based on water systems
        correctionFactor *= (1.0 + 0.65 * (localDensity / 12.0));
        
        // Hydrogen bond correction: enhanced polarization in H-bond networks
        if (inHBondNetwork) {
            correctionFactor *= 1.15;
        }
        
        // Charge asymmetry correction: asymmetric environments need more correction
        correctionFactor *= (1.0 + 0.2 * chargeAsymmetry);
        
        // Apply correction to r0
        r0[i].x *= correctionFactor;
        r0[i].y *= correctionFactor;
        r0[i].z *= correctionFactor;
    }
}

double DrudeForce::calculateLocalDensity(const model::MCState& state, 
                                        const model::MCAtom& drude,
                                        int drudeIdx) {
    double count = 0.0;
    const double cutoff = 0.6;  // 6 Angstrom in nm
    const double cutoff2 = cutoff * cutoff;
    
    // Count nearby oxygen atoms (proxy for water molecules)
    for (int j = 0; j < state.activeAtomCount; j++) {
        const auto& atom = state.atoms[j];
        
        // Look for oxygen atoms (type 0 in our water model)
        if (atom.type != 0) continue;
        
        double dx = drude.x - atom.x;
        double dy = drude.y - atom.y;
        double dz = drude.z - atom.z;
        
        // Apply PBC
        applyPBC(dx, dy, dz, state.info.box[0], state.info.box[1], state.info.box[2]);
        
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < cutoff2 && r2 > 0.01) {
            // Smooth counting function
            double r = std::sqrt(r2);
            count += 0.5 * (1.0 + std::cos(M_PI * r / cutoff));
        }
    }
    
    return count;
}

bool DrudeForce::isInHydrogenBondNetwork(const model::MCState& state,
                                         const DrudeParticle& particle) {
    // Simple H-bond detection based on geometry
    const auto& oxygen = state.atoms[particle.parentIndex];
    
    int hbondCount = 0;
    const double hbondCutoff = 0.35;  // 3.5 Angstrom
    const double hbondCutoff2 = hbondCutoff * hbondCutoff;
    const double angleThreshold = std::cos(30.0 * M_PI / 180.0);  // 30 degrees
    
    // Check all other water molecules
    for (int i = 0; i < state.activeResidueCount; i++) {
        const auto& res = state.residues[i];
        if (res.atomCount != 5) continue;  // Not a water
        
        // Get the other oxygen
        const auto& otherO = state.atoms[res.atomStart];
        if (&otherO == &oxygen) continue;
        
        double dx = oxygen.x - otherO.x;
        double dy = oxygen.y - otherO.y;
        double dz = oxygen.z - otherO.z;
        
        applyPBC(dx, dy, dz, state.info.box[0], state.info.box[1], state.info.box[2]);
        
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < hbondCutoff2 && r2 > 0.01) {
            // Check angle criterion (simplified)
            hbondCount++;
        }
    }
    
    // In H-bond network if participating in 2 or more H-bonds
    return hbondCount >= 2;
}

double DrudeForce::calculateLocalChargeAsymmetry(const model::MCState& state,
                                                 const model::MCAtom& drude,
                                                 int drudeIdx) {
    // Calculate charge distribution asymmetry in local environment
    Vec3 chargeVector(0.0, 0.0, 0.0);
    double totalCharge = 0.0;
    
    const double cutoff = 0.5;  // 5 Angstrom
    const double cutoff2 = cutoff * cutoff;
    
    for (int j = 0; j < state.activeAtomCount; j++) {
        const auto& atom = state.atoms[j];
        if (j == particles[drudeIdx].drudeIndex) continue;  // Skip self
        
        double dx = atom.x - drude.x;
        double dy = atom.y - drude.y;
        double dz = atom.z - drude.z;
        
        applyPBC(dx, dy, dz, state.info.box[0], state.info.box[1], state.info.box[2]);
        
        double r2 = dx*dx + dy*dy + dz*dz;
        if (r2 < cutoff2 && r2 > 0.01) {
            double r = std::sqrt(r2);
            double weight = atom.charge / r2;
            
            chargeVector.x += weight * dx / r;
            chargeVector.y += weight * dy / r;
            chargeVector.z += weight * dz / r;
            
            totalCharge += std::abs(atom.charge) / r;
        }
    }
    
    // Asymmetry metric: magnitude of charge vector normalized by total charge
    double asymmetry = 0.0;
    if (totalCharge > 0) {
        asymmetry = chargeVector.norm() / totalCharge;
    }
    
    return std::min(asymmetry, 1.0);  // Cap at 1.0
}

void DrudeForce::calculateStandardR0(model::MCState& state, std::vector<Vec3>& r0) {
    // Standard r0 calculation (unchanged from original)
    r0.resize(particles.size());
    
    // Reset Drude positions to parent positions
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& particle = particles[i];
        state.atoms[particle.drudeIndex].x = state.atoms[particle.parentIndex].x;
        state.atoms[particle.drudeIndex].y = state.atoms[particle.parentIndex].y;
        state.atoms[particle.drudeIndex].z = state.atoms[particle.parentIndex].z;
    }
    
    // Calculate electric field at each Drude from fixed charges only
    std::vector<Vec3> electricField(particles.size(), Vec3(0.0, 0.0, 0.0));
    calculateElectricFieldAtDrudes(state, electricField, false);  // false = exclude Drudes
    
    // Convert field to displacement
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& particle = particles[i];
        double k = 4184.0 * 1000.0;  // Force constant
        
        r0[i] = electricField[i] * (particle.charge * 138.935 / k);
    }
}

// Note: calculateEnergySmartOPT3 functionality is now handled through 
// minimizeDrudePositionsWithSmartOPT3 which is called by calculateEnergySCF

bool DrudeForce::minimizeDrudePositionsWithSmartOPT3(model::MCState& state) {
    // Calculate smart r0 with corrections
    std::vector<Vec3> r0;
    calculateSmartR0(state, r0);
    
    // Rest follows standard OPT3 algorithm
    std::vector<Vec3> r1, r2, r3;
    
    // Calculate r1
    updateDrudePositions(state, r0);
    calculateElectricFieldResponse(state, r1);
    
    // Calculate r2
    updateDrudePositions(state, r1);
    calculateElectricFieldResponse(state, r2);
    
    // Calculate r3
    updateDrudePositions(state, r2);
    calculateElectricFieldResponse(state, r3);
    
    // Combine with optimized coefficients for smart r0
    // These coefficients are tuned for the corrected r0
    const double c0 = 0.15;   // Now r0 is useful!
    const double c1 = 0.35;
    const double c2 = 0.30;
    const double c3 = 0.20;
    
    // Final positions
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& particle = particles[i];
        Vec3 finalDisp = r0[i] * c0 + r1[i] * c1 + r2[i] * c2 + r3[i] * c3;
        
        state.atoms[particle.drudeIndex].x = state.atoms[particle.parentIndex].x + finalDisp.x;
        state.atoms[particle.drudeIndex].y = state.atoms[particle.parentIndex].y + finalDisp.y;
        state.atoms[particle.drudeIndex].z = state.atoms[particle.parentIndex].z + finalDisp.z;
    }
    
    return true;  // Smart OPT3 doesn't iterate, so always "converges"
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc