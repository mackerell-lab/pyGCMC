#include "EnergyDirectCore.hpp"
#include "../coulomb/CoulombPairCore.hpp"
#include "EnergyUtils.hpp"
#include "platform/Platform.hpp"
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

    const bool includePairtypes14Intra = movement_only;

    if (movement_only) {
        for (const auto& movementInfo : state.movementResidues) {
            for (int i = movementInfo.startIndex;
                 i < movementInfo.startIndex + movementInfo.activeCount;
                 ++i) {
                if (!residues[i].active) continue;
                computeResidueNonbondedEnergy(state, i, use_cutoff, use_pbc,
                                              vdw_only, includePairtypes14Intra);
            }
        }
    } else {
        for (int i = 0; i < state.activeResidueCount; ++i) {
            if (!residues[i].active) continue;
            computeResidueNonbondedEnergy(state, i, use_cutoff, use_pbc,
                                          vdw_only, includePairtypes14Intra);
        }
    }
}

void computeResidueNonbondedEnergy(model::MCState& state,
                                   int residue_idx,
                                   bool use_cutoff,
                                   bool use_pbc,
                                   bool vdw_only,
                                   bool include_pairtypes14_intra,
                                   ResiduePartnerFilter partner_filter) {
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
    const bool usePairtypes14 = forcefield.pairtypes14Enabled;

    if (!residues[residue_idx].active) {
        return;
    }

    residues[residue_idx].energy_vdw = 0.0f;
    residues[residue_idx].energy_elec = 0.0f;

    for (int atom_i = residues[residue_idx].atomStart;
         atom_i < residues[residue_idx].atomStart + residues[residue_idx].atomCount;
         ++atom_i) {
        // Skip LP/LPA atoms (lone pairs)
        if (atoms[atom_i].name == "LP" || atoms[atom_i].name == "LPA") continue;

        int type_i = atoms[atom_i].type;

        for (int j = 0; j < state.activeResidueCount; ++j) {
            if (!residues[j].active || j == residue_idx) continue;
            if (partner_filter == ResiduePartnerFilter::FixedOnly && !residues[j].fixed) continue;
            if (partner_filter == ResiduePartnerFilter::NonFixedOnly && residues[j].fixed) continue;

            for (int atom_j = residues[j].atomStart;
                 atom_j < residues[j].atomStart + residues[j].atomCount;
                 ++atom_j) {
                // Skip LP/LPA atoms (lone pairs)
                if (atoms[atom_j].name == "LP" || atoms[atom_j].name == "LPA") continue;

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

    if (include_pairtypes14_intra && usePairtypes14) {
        const int start = residues[residue_idx].atomStart;
        const int end = start + residues[residue_idx].atomCount;
        const int nTypes = forcefield.numTotalTypes;

        for (int atom_i = start; atom_i < end; ++atom_i) {
            if (atoms[atom_i].name == "LP" || atoms[atom_i].name == "LPA") continue;
            const int type_i = atoms[atom_i].type;

            for (int atom_j = atom_i + 1; atom_j < end; ++atom_j) {
                if (atoms[atom_j].name == "LP" || atoms[atom_j].name == "LPA") continue;
                if (!state.isPair14(atom_i, atom_j)) continue;

                const int type_j = atoms[atom_j].type;
                if (!forcefield.hasPairtype14(type_i, type_j)) continue;

                double dx = atoms[atom_j].x - atoms[atom_i].x;
                double dy = atoms[atom_j].y - atoms[atom_i].y;
                double dz = atoms[atom_j].z - atoms[atom_i].z;

                if (use_pbc) {
                    dx -= box[0] * std::round(dx / box[0]);
                    dy -= box[1] * std::round(dy / box[1]);
                    dz -= box[2] * std::round(dz / box[2]);
                }

                const double r2 = dx * dx + dy * dy + dz * dz;
                if (r2 > cutoff2 || r2 <= 0.0) continue;

                const int param_index = type_i * nTypes + type_j;
                const size_t idx = static_cast<size_t>(param_index);
                if (idx >= forcefield.ljSigma14.size() || idx >= forcefield.ljEps14.size()) {
                    continue;
                }

                const double sigma = forcefield.ljSigma14[idx];
                const double eps = forcefield.ljEps14[idx];
                const double inv_r2 = 1.0 / r2;
                const double sr2 = (sigma * sigma) * inv_r2;
                const double sr6 = sr2 * sr2 * sr2;
                const double sr12 = sr6 * sr6;
                const double vdw = 4.0 * eps * (sr12 - sr6);
                residues[residue_idx].energy_vdw += static_cast<float>(vdw);
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
    computeResidueNonbondedEnergy(state, residue_idx, true, true, false, true);
}

void computeNonbondedEnergyWithNeighborList(
    model::MCState& state,
    const NeighborList& neighborList,
    bool use_pbc,
    bool vdw_only
) {
    if (neighborList.empty()) {
        // Fallback to full calculation if neighbor list not available
        computeNonbondedEnergy(state, true, false, use_pbc, vdw_only);
        return;
    }

    auto& residues = state.residues;
    auto& forcefield = state.forcefield;

    // Build NxN matrix if needed
    const size_t expected_size_check = static_cast<size_t>(forcefield.numTotalTypes) *
                                       static_cast<size_t>(forcefield.numTotalTypes);

    if ((!forcefield.ljMatrixInitialized) ||
        (forcefield.ljSigma.size() != expected_size_check) ||
        (forcefield.ljEps.size() != expected_size_check)) {
        forcefield.rebuildLJMatrix();
    }

    const auto& atoms = state.atoms;
    const auto& box = state.info.box;

    // Build atom to residue mapping
    std::vector<int> atom_to_residue(state.activeAtomCount, -1);
    for (int res_idx = 0; res_idx < state.activeResidueCount; ++res_idx) {
        if (!residues[res_idx].active) continue;
        for (int atom_idx = residues[res_idx].atomStart;
             atom_idx < residues[res_idx].atomStart + residues[res_idx].atomCount;
             ++atom_idx) {
            if (atom_idx >= 0 && atom_idx < state.activeAtomCount) {
                atom_to_residue[atom_idx] = res_idx;
            }
        }
    }

    // Reset energies for all residues
    for (auto& residue : residues) {
        residue.energy_vdw = 0.0f;
        if (!vdw_only) {
            residue.energy_elec = 0.0f;
        }
    }

    // Use neighbor list to compute pair energies
    for (int atom_i = 0; atom_i < state.activeAtomCount; ++atom_i) {
        // Skip LP/LPA atoms (lone pairs)
        if (atoms[atom_i].name == "LP" || atoms[atom_i].name == "LPA") continue;

        int type_i = atoms[atom_i].type;
        int res_i = atom_to_residue[atom_i];

        if (res_i < 0 || res_i >= static_cast<int>(residues.size()) || !residues[res_i].active) {
            continue;
        }

        const auto& neighbors = neighborList.getNeighbors(atom_i);

        for (int atom_j : neighbors) {
            // Only compute each pair once (atom_j > atom_i ensures this)
            if (atom_j <= atom_i || atom_j >= state.activeAtomCount) continue;

            // Skip LP/LPA atoms
            if (atoms[atom_j].name == "LP" || atoms[atom_j].name == "LPA") continue;

            int type_j = atoms[atom_j].type;
            int res_j = atom_to_residue[atom_j];

            if (res_j < 0 || res_j >= static_cast<int>(residues.size()) || !residues[res_j].active) {
                continue;
            }

            // Skip intra-residue interactions
            if (res_i == res_j) continue;

            double dx = atoms[atom_j].x - atoms[atom_i].x;
            double dy = atoms[atom_j].y - atoms[atom_i].y;
            double dz = atoms[atom_j].z - atoms[atom_i].z;

            if (use_pbc) {
                dx -= box[0] * std::round(dx / box[0]);
                dy -= box[1] * std::round(dy / box[1]);
                dz -= box[2] * std::round(dz / box[2]);
            }

            double r2 = dx*dx + dy*dy + dz*dz;

            int param_index = type_i * forcefield.numTotalTypes + type_j;
            size_t idx = static_cast<size_t>(param_index);

            if (idx >= forcefield.ljEps.size() || idx >= forcefield.ljSigma.size()) {
                continue;
            }

            double eps = forcefield.ljEps[idx];
            double sigma = forcefield.ljSigma[idx];
            double q1 = atoms[atom_i].charge;
            double q2 = atoms[atom_j].charge;

            auto [vdw, elec] = coulomb::calcPairEnergy(r2, sigma, eps, q1, q2, state.info, !vdw_only);

            // Add energy to both residues (each pair counted once)
            residues[res_i].energy_vdw += static_cast<float>(vdw);
            residues[res_j].energy_vdw += static_cast<float>(vdw);
            if (!vdw_only) {
                residues[res_i].energy_elec += static_cast<float>(elec);
                residues[res_j].energy_elec += static_cast<float>(elec);
            }
        }
    }
}

void computeResidueEnergyWithNeighborList(
    model::MCState& state,
    int residue_idx,
    const NeighborList& neighborList,
    bool use_pbc,
    bool vdw_only
) {
    if (neighborList.empty()) {
        // Fallback to full calculation
        computeResidueNonbondedEnergy(state, residue_idx, true, use_pbc, vdw_only, false);
        return;
    }

    auto& residues = state.residues;
    auto& forcefield = state.forcefield;

    // Build NxN matrix if needed
    const size_t expected_size_check = static_cast<size_t>(forcefield.numTotalTypes) *
                                       static_cast<size_t>(forcefield.numTotalTypes);

    if ((!forcefield.ljMatrixInitialized) ||
        (forcefield.ljSigma.size() != expected_size_check) ||
        (forcefield.ljEps.size() != expected_size_check)) {
        forcefield.rebuildLJMatrix();
    }

    const auto& atoms = state.atoms;
    const auto& box = state.info.box;

    if (!residues[residue_idx].active) {
        return;
    }

    residues[residue_idx].energy_vdw = 0.0f;
    residues[residue_idx].energy_elec = 0.0f;

    for (int atom_i = residues[residue_idx].atomStart;
         atom_i < residues[residue_idx].atomStart + residues[residue_idx].atomCount;
         ++atom_i) {
        // Skip LP/LPA atoms
        if (atoms[atom_i].name == "LP" || atoms[atom_i].name == "LPA") continue;

        int type_i = atoms[atom_i].type;

        const auto& neighbors = neighborList.getNeighbors(atom_i);

        for (int atom_j : neighbors) {
            if (atom_j < 0 || atom_j >= state.activeAtomCount) continue;

            // Skip LP/LPA atoms
            if (atoms[atom_j].name == "LP" || atoms[atom_j].name == "LPA") continue;

            // Find which residue atom_j belongs to
            int res_j = -1;
            for (int r = 0; r < state.activeResidueCount; ++r) {
                if (!residues[r].active) continue;
                if (atom_j >= residues[r].atomStart &&
                    atom_j < residues[r].atomStart + residues[r].atomCount) {
                    res_j = r;
                    break;
                }
            }

            if (res_j < 0 || res_j >= static_cast<int>(residues.size()) || !residues[res_j].active) {
                continue;
            }

            // Skip intra-residue interactions
            if (res_j == residue_idx) continue;

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

            int param_index = type_i * forcefield.numTotalTypes + type_j;
            size_t idx = static_cast<size_t>(param_index);

            if (idx >= forcefield.ljEps.size() || idx >= forcefield.ljSigma.size()) {
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

} // namespace cpu
} // namespace platform
} // namespace pygcmc
