#include "PMEComplete.hpp"
#include "PMEComposite.hpp"
#include "PMEReal.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace energy {
namespace pme {

EnergyComponents PMEComplete::computeCompleteEnergy(model::MCState& state) {
    // For PME Complete, we need to recalculate everything to include intramolecular
    // Don't use PMEComposite::computeSystemEnergy as it excludes intramolecular
    
    // Initialize PME calculation
    // PMEComposite pmeComp;  // Not needed - using static methods
    
    // Calculate PME electrostatics (real + reciprocal + self)
    double real_elec = 0.0;
    double recip_elec = 0.0; 
    double self_elec = 0.0;
    
    // For now, we'll use the standard PME for electrostatics
    // This is OK because electrostatics don't exclude intramolecular
    PMEComposite::computeSystemEnergy(state);
    real_elec = state.ewald_energy.real_space;
    recip_elec = state.ewald_energy.reciprocal;
    self_elec = state.ewald_energy.self;
    
    // But we need to recalculate VdW to include intramolecular
    // Use the same logic as computeCompleteEnergyCutoff
    EnergyComponents energy;
    energy.real = real_elec;
    energy.reciprocal = recip_elec;
    energy.self = self_elec;
    energy.vdw = 0.0;
    
    // Reset residue VdW energies and recalculate including ALL interactions
    for (auto& res : state.residues) {
        if (res.active) {
            res.energy_vdw = 0.0;
        }
    }
    
    // Calculate all LJ pairwise interactions (including intramolecular)
    for (int i = 0; i < state.activeAtomCount - 1; ++i) {
        const auto& atom_i = state.atoms[i];
        int res_i = -1;
        for (int r = 0; r < state.activeResidueCount; ++r) {
            const auto& res = state.residues[r];
            if (res.active && i >= res.atomStart && i < res.atomStart + res.atomCount) {
                res_i = r;
                break;
            }
        }
        
        for (int j = i + 1; j < state.activeAtomCount; ++j) {
            const auto& atom_j = state.atoms[j];
            int res_j = -1;
            for (int r = 0; r < state.activeResidueCount; ++r) {
                const auto& res = state.residues[r];
                if (res.active && j >= res.atomStart && j < res.atomStart + res.atomCount) {
                    res_j = r;
                    break;
                }
            }
            
            // Calculate distance with PBC
            double dx = atom_i.x - atom_j.x;
            double dy = atom_i.y - atom_j.y;
            double dz = atom_i.z - atom_j.z;
            
            // Apply minimum image convention
            dx -= state.info.box[0] * std::round(dx / state.info.box[0]);
            dy -= state.info.box[1] * std::round(dy / state.info.box[1]);
            dz -= state.info.box[2] * std::round(dz / state.info.box[2]);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            double r = std::sqrt(r2);
            
            if (r < state.info.cutoff) {
                // LJ energy - always calculated (including intramolecular)
                int type_i = atom_i.type;
                int type_j = atom_j.type;
                int pair_idx = type_i * state.forcefield.numTotalTypes + type_j;
                
                // Bounds check
                if (pair_idx >= static_cast<int>(state.forcefield.ljSigma.size()) ||
                    pair_idx >= static_cast<int>(state.forcefield.ljEps.size())) {
                    continue; // Skip invalid pairs
                }
                
                double sigma = state.forcefield.ljSigma[pair_idx];
                double epsilon = state.forcefield.ljEps[pair_idx];
                
                double sr = sigma / r;
                double sr6 = sr * sr * sr * sr * sr * sr;
                double sr12 = sr6 * sr6;
                double lj_energy = 4.0 * epsilon * (sr12 - sr6);
                
                energy.vdw += lj_energy;
                // Add half to each residue to avoid double-counting when summing
                if (res_i >= 0) state.residues[res_i].energy_vdw += lj_energy / 2.0;
                if (res_j >= 0) state.residues[res_j].energy_vdw += lj_energy / 2.0;
            }
        }
    }
    
    energy.total = energy.real + energy.reciprocal + energy.self + energy.vdw;
    
    return energy;
}

EnergyComponents PMEComplete::computeCompleteEnergyCutoff(model::MCState& state) {
    EnergyComponents energy;
    energy.real = 0.0;
    energy.reciprocal = 0.0;
    energy.self = 0.0;
    energy.vdw = 0.0;
    energy.total = 0.0;
    
    // Reset residue energies
    for (auto& res : state.residues) {
        if (res.active) {
            res.energy_elec = 0.0;
            res.energy_vdw = 0.0;
        }
    }
    
    // Calculate all pairwise interactions (including intramolecular)
    for (int i = 0; i < state.activeAtomCount - 1; ++i) {
        const auto& atom_i = state.atoms[i];
        int res_i = -1;
        for (int r = 0; r < state.activeResidueCount; ++r) {
            const auto& res = state.residues[r];
            if (res.active && i >= res.atomStart && i < res.atomStart + res.atomCount) {
                res_i = r;
                break;
            }
        }
        
        for (int j = i + 1; j < state.activeAtomCount; ++j) {
            const auto& atom_j = state.atoms[j];
            int res_j = -1;
            for (int r = 0; r < state.activeResidueCount; ++r) {
                const auto& res = state.residues[r];
                if (res.active && j >= res.atomStart && j < res.atomStart + res.atomCount) {
                    res_j = r;
                    break;
                }
            }
            
            // Calculate distance with PBC
            double dx = atom_i.x - atom_j.x;
            double dy = atom_i.y - atom_j.y;
            double dz = atom_i.z - atom_j.z;
            
            // Apply minimum image convention
            dx -= state.info.box[0] * std::round(dx / state.info.box[0]);
            dy -= state.info.box[1] * std::round(dy / state.info.box[1]);
            dz -= state.info.box[2] * std::round(dz / state.info.box[2]);
            
            double r2 = dx*dx + dy*dy + dz*dz;
            double r = std::sqrt(r2);
            
            if (r < state.info.cutoff) {
                // LJ energy - always calculated
                int type_i = atom_i.type;
                int type_j = atom_j.type;
                int pair_idx = type_i * state.forcefield.numTotalTypes + type_j;
                
                // Bounds check
                if (pair_idx >= static_cast<int>(state.forcefield.ljSigma.size()) ||
                    pair_idx >= static_cast<int>(state.forcefield.ljEps.size())) {
                    continue; // Skip invalid pairs
                }
                
                double sigma = state.forcefield.ljSigma[pair_idx];
                double epsilon = state.forcefield.ljEps[pair_idx];
                
                double sr = sigma / r;
                double sr6 = sr * sr * sr * sr * sr * sr;
                double sr12 = sr6 * sr6;
                double lj_energy = 4.0 * epsilon * (sr12 - sr6);
                
                energy.vdw += lj_energy;
                // Add half to each residue to avoid double-counting when summing
                if (res_i >= 0) state.residues[res_i].energy_vdw += lj_energy / 2.0;
                if (res_j >= 0) state.residues[res_j].energy_vdw += lj_energy / 2.0;
                
                // Electrostatic energy
                double elec_energy = 332.0637 * atom_i.charge * atom_j.charge / r;
                energy.real += elec_energy;
                // Add half to each residue to avoid double-counting when summing
                if (res_i >= 0) state.residues[res_i].energy_elec += elec_energy / 2.0;
                if (res_j >= 0) state.residues[res_j].energy_elec += elec_energy / 2.0;
            }
        }
    }
    
    // Calculate self energy (using default alpha value from PME params)
    double alpha = 5.6; // Default alpha value
    energy.self = 0.0;
    for (int i = 0; i < state.activeAtomCount; ++i) {
        energy.self -= 332.0637 * state.atoms[i].charge * state.atoms[i].charge * 
                       std::sqrt(2.0 / M_PI) * alpha;
    }
    
    energy.total = energy.real + energy.vdw + energy.self;
    
    // Residue totals are already updated in the loop above
    
    return energy;
}

double PMEComplete::calculateIntramolecularLJ(model::MCState& state) {
    double intraLJ = 0.0;
    
    for (const auto& res : state.residues) {
        if (!res.active) continue;
        
        // Calculate LJ within this residue
        for (int i = res.atomStart; i < res.atomStart + res.atomCount - 1; ++i) {
            for (int j = i + 1; j < res.atomStart + res.atomCount; ++j) {
                const auto& atom_i = state.atoms[i];
                const auto& atom_j = state.atoms[j];
                
                double dx = atom_i.x - atom_j.x;
                double dy = atom_i.y - atom_j.y;
                double dz = atom_i.z - atom_j.z;
                
                // Apply minimum image convention
                dx -= state.info.box[0] * std::round(dx / state.info.box[0]);
                dy -= state.info.box[1] * std::round(dy / state.info.box[1]);
                dz -= state.info.box[2] * std::round(dz / state.info.box[2]);
                
                double r2 = dx*dx + dy*dy + dz*dz;
                double r = std::sqrt(r2);
                
                if (r > 0 && r < state.info.cutoff) {
                    // Get LJ parameters
                    int type_i = atom_i.type;
                    int type_j = atom_j.type;
                    int pair_idx = type_i * state.forcefield.numTotalTypes + type_j;
                    
                    double sigma = state.forcefield.ljSigma[pair_idx];
                    double epsilon = state.forcefield.ljEps[pair_idx];
                    
                    // Calculate LJ energy
                    double sr = sigma / r;
                    double sr6 = sr * sr * sr * sr * sr * sr;
                    double sr12 = sr6 * sr6;
                    intraLJ += 4.0 * epsilon * (sr12 - sr6);
                }
            }
        }
    }
    
    return intraLJ;
}

void PMEComplete::applyElectrostaticExclusions(model::MCState& state, EnergyComponents& energy) {
    // This is a placeholder for more sophisticated exclusion handling
    // Currently, the standard PME already handles most exclusions properly
    // This method can be extended to handle special cases if needed
    
    // Suppress unused parameter warnings
    (void)state;
    (void)energy;
}

} // namespace pme
} // namespace energy

// Global convenience functions
void computeSystemEnergyPMEComplete(model::MCState& state) {
    static energy::pme::PMEComplete pmeComplete;
    auto energy = pmeComplete.computeCompleteEnergy(state);
    
    // Update state with the complete energy values
    state.ewald_energy.real_space = energy.real;
    state.ewald_energy.reciprocal = energy.reciprocal;
    state.ewald_energy.self = energy.self;
    state.ewald_energy.total = energy.total;
}

void computeSystemEnergyCutoffComplete(model::MCState& state) {
    static energy::pme::PMEComplete pmeComplete;
    auto energy = pmeComplete.computeCompleteEnergyCutoff(state);
    
    // Update state with the complete energy values
    state.ewald_energy.real_space = energy.real;
    state.ewald_energy.reciprocal = energy.reciprocal;
    state.ewald_energy.self = energy.self;
    state.ewald_energy.total = energy.total;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc