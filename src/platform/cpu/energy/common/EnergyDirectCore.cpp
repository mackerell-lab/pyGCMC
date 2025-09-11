#include "EnergyDirectCore.hpp"
#include "../coulomb/CoulombPairCore.hpp"
#include "EnergyUtils.hpp"
#include "platform/platform.hpp"
#include <stdexcept>
#include <sstream>
#include <limits>
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

void computeNonbondedEnergy(model::MCState& state, bool use_cutoff, bool movement_only, bool use_pbc, bool vdw_only) {
    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << "\n=== Starting nonbonded energy calculation ===";
        ss << "\nSystem state info:";
        ss << "\n  Active residue count: " << state.activeResidueCount;
        if (use_cutoff) {
            ss << "\n  Cutoff distance: " << state.info.cutoff << " nm";
        }
        if (movement_only) {
            ss << "\n  Calculating only for movement residues";
        }
        if (use_pbc) {
            ss << "\n  Using periodic boundary conditions";
        }
        if (vdw_only) {
            ss << "\n  Calculating VDW interactions only";
        }
        platform::log(LogLevel::DEBUG, ss.str());
    }

    auto& residues = state.residues;
    auto& forcefield = state.forcefield;  // Non-const to allow rebuild
    
    // Build NxN matrix if needed (from mixing rules + NBFIX)
    const size_t expected_size_check = static_cast<size_t>(forcefield.numTotalTypes) * 
                                       static_cast<size_t>(forcefield.numTotalTypes);
    
    if ((!forcefield.ljMatrixInitialized) ||
        (forcefield.ljSigma.size() != expected_size_check) ||
        (forcefield.ljEps.size() != expected_size_check)) {
        forcefield.rebuildLJMatrix();
    }

    // Validate basic state parameters
    if (state.activeResidueCount < 0 || 
        static_cast<size_t>(state.activeResidueCount) > residues.size()) {
        throw std::runtime_error("Invalid activeResidueCount: " + 
                               std::to_string(state.activeResidueCount) +
                               " (residues size: " + std::to_string(residues.size()) + ")");
    }

    if (forcefield.numTotalTypes <= 0) {
        throw std::runtime_error("Invalid numTotalTypes: " + 
                               std::to_string(forcefield.numTotalTypes));
    }

    if (movement_only && forcefield.numMovementTypes <= 0) {
        throw std::runtime_error("Invalid numMovementTypes: " + 
                               std::to_string(forcefield.numMovementTypes));
    }

    // Always expect full matrix size
    size_t expected_size = static_cast<size_t>(forcefield.numTotalTypes) * 
                           static_cast<size_t>(forcefield.numTotalTypes);
    
    if (forcefield.ljEps.size() != expected_size || forcefield.ljSigma.size() != expected_size) {
        std::stringstream ss;
        ss << "Force field parameters array size mismatch. Expected size "
           << expected_size
           << " (numTotalTypes * numTotalTypes), but got eps=" << forcefield.ljEps.size()
           << " sigma=" << forcefield.ljSigma.size();
        throw std::runtime_error(ss.str());
    }

    // Reset energies for all residues
    for (auto& residue : residues) {
        residue.energy_vdw = 0.0f;
        if (!vdw_only) {
            residue.energy_elec = 0.0f;
        }
    }

    if (movement_only) {
        for (const auto& movementInfo : state.movementResidues) {
            for (int i = movementInfo.startIndex;
                 i < movementInfo.startIndex + movementInfo.activeCount;
                 ++i) {
                if (!residues[i].active) continue;
                computeResidueNonbondedEnergy(state, i, use_cutoff, use_pbc, vdw_only);
            }
        }
    } else {
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (!residues[i].active) continue;
            computeResidueNonbondedEnergy(state, i, use_cutoff, use_pbc, vdw_only);
        }
    }
}

void computeResidueNonbondedEnergy(model::MCState& state, int residue_idx, bool use_cutoff, bool use_pbc, bool vdw_only) {
    auto& residues = state.residues;
    auto& forcefield = state.forcefield;  // Non-const to allow rebuild
    
    // Build NxN matrix if needed (from mixing rules + NBFIX)
    const size_t expected_size_check = static_cast<size_t>(forcefield.numTotalTypes) * 
                                       static_cast<size_t>(forcefield.numTotalTypes);
    
    if ((!forcefield.ljMatrixInitialized) ||
        (forcefield.ljSigma.size() != expected_size_check) ||
        (forcefield.ljEps.size() != expected_size_check)) {
        forcefield.rebuildLJMatrix();
    }
    const auto& atoms = state.atoms;
    const auto& box = state.info.box;
    
    const double cutoff2 = use_cutoff ? state.info.cutoff * state.info.cutoff : std::numeric_limits<double>::max();
    
    if (!residues[residue_idx].active) {
        return;
    }
    
    residues[residue_idx].energy_vdw = 0.0f;
    residues[residue_idx].energy_elec = 0.0f;
    
    for (int atom_i = residues[residue_idx].atomStart;
         atom_i < residues[residue_idx].atomStart + residues[residue_idx].atomCount;
         ++atom_i) {
        int type_i = atoms[atom_i].type;
        
        for (int j = 0; j < state.activeResidueCount; ++j) {
            if (!residues[j].active || j == residue_idx) continue;
            
            for (int atom_j = residues[j].atomStart;
                 atom_j < residues[j].atomStart + residues[j].atomCount;
                 ++atom_j) {
                int type_j = atoms[atom_j].type;
                
                double dx = atoms[atom_j].x - atoms[atom_i].x;
                double dy = atoms[atom_j].y - atoms[atom_i].y;
                double dz = atoms[atom_j].z - atoms[atom_i].z;
                
                if (use_pbc) {
                    dx -= box[0] * std::round(dx / box[0]);
                    dy -= box[1] * std::round(dy / box[1]);
                    dz -= box[2] * std::round(dz / box[2]);
                }
                
                double r2 = dx*dx + dy*dy + dz*dz;
                
                if (r2 > cutoff2) continue;
                
                int param_index = type_i * forcefield.numTotalTypes + type_j;
                size_t idx = static_cast<size_t>(param_index);
                
                // Bounds check to prevent accessing invalid force field parameters
                if (idx >= forcefield.ljEps.size() || idx >= forcefield.ljSigma.size()) {
                    // Skip if force field tables don't cover this pair
                    continue;
                }
                
                double eps = forcefield.ljEps[idx];
                double sigma = forcefield.ljSigma[idx];
                double q1 = atoms[atom_i].charge;
                double q2 = atoms[atom_j].charge;
                
                auto [vdw, elec] = coulomb::calcPairEnergy(r2, sigma, eps, q1, q2, state.info, !vdw_only);
                
                residues[residue_idx].energy_vdw += static_cast<float>(vdw);
                if (!vdw_only) {
                    residues[residue_idx].energy_elec += static_cast<float>(elec);
                }
            }
        }
    }
}

// Basic functions
void computeMovementEnergy(model::MCState& state) {
    computeNonbondedEnergy(state, false, true, false);
}

void computeMovementEnergyCutoff(model::MCState& state) {
    computeNonbondedEnergy(state, true, true, false);
}

void computeSystemEnergy(model::MCState& state) {
    computeNonbondedEnergy(state, false, false, false);
}

void computeSystemEnergyCutoff(model::MCState& state) {
    computeNonbondedEnergy(state, true, false, false);
}

void computeSystemVdwEnergyCutoff(model::MCState& state) {
    computeNonbondedEnergy(state, true, false, false, true);
}

void computeSystemEnergyPBC(model::MCState& state) {
    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << "\n=== Starting PBC nonbonded energy calculation (no cutoff) ===";
        ss << "\nBox dimensions: " << state.info.box[0] << " x " 
           << state.info.box[1] << " x " << state.info.box[2] << " nm";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // Validate box dimensions before proceeding
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        throw std::runtime_error("Invalid box dimensions for PBC calculation");
    }
    
    computeNonbondedEnergy(state, false, false, true);
}

void computeSystemEnergyPBCCutoff(model::MCState& state) {
    if (getEnergyDebugOutput()) {
        std::stringstream ss;
        ss << "\n=== Starting PBC nonbonded energy calculation (with cutoff) ===";
        ss << "\nBox dimensions: " << state.info.box[0] << " x " 
           << state.info.box[1] << " x " << state.info.box[2] << " nm";
        platform::log(LogLevel::DEBUG, ss.str());
    }
    
    // Validate box dimensions before proceeding
    if (state.info.box[0] <= 0.0f || state.info.box[1] <= 0.0f || state.info.box[2] <= 0.0f) {
        throw std::runtime_error("Invalid box dimensions for PBC calculation");
    }
    
    computeNonbondedEnergy(state, true, false, true);
}

void computeSystemEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc) {
    computeNonbondedEnergy(state, use_cutoff, false, use_pbc);
}

void computeMovementEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc) {
    computeNonbondedEnergy(state, use_cutoff, true, use_pbc);
}

void computeSystemVdwEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc) {
    computeNonbondedEnergy(state, use_cutoff, false, use_pbc, true);
}

void computeMovementVdwEnergyDirect(model::MCState& state, bool use_cutoff, bool use_pbc) {
    computeNonbondedEnergy(state, use_cutoff, true, use_pbc, true);
}

void computeResidueEnergyCutoffPBC(model::MCState& state, int residue_idx) {
    // Wrapper for multi-insertion optimization: single residue energy with PBC and cutoff
    computeResidueNonbondedEnergy(state, residue_idx, true, true, false);
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 