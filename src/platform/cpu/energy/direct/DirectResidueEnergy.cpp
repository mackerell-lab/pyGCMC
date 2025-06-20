#include "DirectResidueEnergy.hpp"
#include "../common/EnergyUtils.hpp"
#include "platform/platform.hpp"
#include <stdexcept>
#include <sstream>
#include <limits>
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace direct {

/**
 * @brief Calculate nonbonded interactions between a single residue and all other active residues
 */
void computeResidueNonbondedEnergy(
    model::MCState& state,
    int residue_idx,
    bool use_cutoff,
    bool use_pbc
) {
    auto& residues = state.residues;
    const auto& forcefield = state.forcefield;
    const auto& atoms = state.atoms;
    const auto& box = state.info.box;  // Box dimensions for PBC
    
    // Calculate squared cutoff distance if using cutoff
    const double cutoff2 = use_cutoff ? state.info.cutoff * state.info.cutoff : std::numeric_limits<double>::max();
    
    // Return if residue is not active
    if (!residues[residue_idx].active) {
        return;
    }
    
    // Reset energy components for current residue
    residues[residue_idx].energy_vdw = 0.0f;
    residues[residue_idx].energy_elec = 0.0f;
    
    // Iterate through all atoms in current residue
    for (int atom_i = residues[residue_idx].atomStart;
         atom_i < residues[residue_idx].atomStart + residues[residue_idx].atomCount;
         ++atom_i) {
        int type_i = atoms[atom_i].type;
        
        // Validate atom type
        if (type_i >= forcefield.numTotalTypes) {
            std::stringstream ss;
            ss << "Atom type " << type_i << " out of range. "
               << "Maximum allowed type is " << (forcefield.numTotalTypes - 1)
               << " for atom " << atom_i << " in residue " << residue_idx;
            throw std::runtime_error(ss.str());
        }
        
        // Calculate interactions with atoms in other active residues
        for (int j = 0; j < state.activeResidueCount; ++j) {
            if (!residues[j].active || j == residue_idx) continue;
            
            // Iterate through atoms in other residue
            for (int atom_j = residues[j].atomStart;
                 atom_j < residues[j].atomStart + residues[j].atomCount;
                 ++atom_j) {
                int type_j = atoms[atom_j].type;
                
                // Validate atom type
                if (type_j >= forcefield.numTotalTypes) {
                    std::stringstream ss;
                    ss << "Atom type " << type_j << " out of range. "
                       << "Maximum allowed type is " << (forcefield.numTotalTypes - 1)
                       << " for atom " << atom_j << " in residue " << j;
                    throw std::runtime_error(ss.str());
                }
                
                // Calculate interatomic distance with PBC if enabled
                double dx = atoms[atom_j].x - atoms[atom_i].x;
                double dy = atoms[atom_j].y - atoms[atom_i].y;
                double dz = atoms[atom_j].z - atoms[atom_i].z;
                
                // Apply minimum image convention if PBC is enabled
                if (use_pbc) {
                    // Validate box dimensions
                    if (box[0] <= 0.0f || box[1] <= 0.0f || box[2] <= 0.0f) {
                        throw std::runtime_error("Invalid box dimensions for PBC calculation");
                    }
                    
                    // Apply minimum image convention
                    dx -= box[0] * std::round(dx / box[0]);
                    dy -= box[1] * std::round(dy / box[1]);
                    dz -= box[2] * std::round(dz / box[2]);
                    
                    if (getEnergyDebugOutput()) {
                        std::stringstream ss;
                        ss << "\nPBC distance calculation:";
                        ss << "\n  Original dx,dy,dz: " << (atoms[atom_j].x - atoms[atom_i].x)
                           << ", " << (atoms[atom_j].y - atoms[atom_i].y)
                           << ", " << (atoms[atom_j].z - atoms[atom_i].z);
                        ss << "\n  After PBC dx,dy,dz: " << dx << ", " << dy << ", " << dz;
                        ss << "\n  Box dimensions: " << box[0] << ", " << box[1] << ", " << box[2];
                        platform::log(LogLevel::DEBUG, ss.str());
                    }
                }
                
                double r2 = dx*dx + dy*dy + dz*dz;
                
                // Skip if beyond cutoff distance
                if (r2 > cutoff2) continue;
                
                // Get force field parameters
                int param_index = type_i * forcefield.numTotalTypes + type_j;
                double eps = forcefield.ljEps[param_index];
                double sigma = forcefield.ljSigma[param_index];
                double q1 = atoms[atom_i].charge;
                double q2 = atoms[atom_j].charge;
                
                // Calculate pair energy
                auto [vdw, elec] = calcPairEnergy(r2, sigma, eps, q1, q2, state.info, true);
                
                // Add energy components to current residue
                residues[residue_idx].energy_vdw += static_cast<float>(vdw);
                residues[residue_idx].energy_elec += static_cast<float>(elec);
            }
        }
    }
}

} // namespace direct
} // namespace cpu
} // namespace platform
} // namespace pygcmc 