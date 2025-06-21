#include "MCCore.hpp"
#include <algorithm>
#include <stdexcept>
#include <cmath>

namespace pygcmc {
namespace system {
namespace montecarlo {

void MCCore::addInitialResidues(model::MCState& state, const model::MCResidue* resVec, int resCount,
                               const model::MCAtom* atomVec, int atomCount) {
    if (resCount > state.info.maxResidues || atomCount > state.info.maxAtoms) {
        throw std::runtime_error("Initial system exceeds max capacity");
    }

    std::copy(resVec, resVec + resCount, state.residues.begin());
    std::copy(atomVec, atomVec + atomCount, state.atoms.begin());

    state.activeResidueCount = resCount;
    state.activeAtomCount = atomCount;
}

int MCCore::insertResidue(model::MCState& state, const model::MCResidue& res, const model::MCAtom* atoms) {
    if (state.activeResidueCount >= state.info.maxResidues ||
        state.activeAtomCount + res.atomCount > state.info.maxAtoms) {
        return -1;
    }

    int resIdx = state.activeResidueCount++;
    state.residues[resIdx] = res;
    state.residues[resIdx].atomStart = state.activeAtomCount;
    state.residues[resIdx].active = true;

    std::copy(atoms, atoms + res.atomCount, state.atoms.begin() + state.activeAtomCount);
    state.activeAtomCount += res.atomCount;

    return resIdx;
}

bool MCCore::removeResidue(model::MCState& state, int resIdx) {
    if (resIdx < 0 || resIdx >= state.activeResidueCount || !state.residues[resIdx].active) {
        return false;
    }

    // Get residue info and mark it as inactive
    model::MCResidue& res = state.residues[resIdx];
    res.active = false;
    int atomStart = res.atomStart;
    int atomCount = res.atomCount;

    // Move atoms if needed
    if (atomStart != state.activeAtomCount - atomCount) {
        for (int i = 0; i < atomCount; ++i) {
            state.atoms[atomStart + i] = state.atoms[state.activeAtomCount - atomCount + i];
        }
    }
    state.activeAtomCount -= atomCount;

    // If this is not the last residue, move the last active one to this position
    if (resIdx != state.activeResidueCount - 1) {
        state.residues[resIdx] = state.residues[state.activeResidueCount - 1];
        // Update the moved residue's atom start position if needed
        if (atomStart != state.activeAtomCount) {
            state.residues[resIdx].atomStart = atomStart;
        }
    }
    
    state.activeResidueCount--;
    return true;
}

void MCCore::translateResidue(model::MCState& state, int resIdx, float dx, float dy, float dz) {
    if (resIdx >= 0 && resIdx < state.activeResidueCount) {
        model::MCResidue& res = state.residues[resIdx];
        // dx, dy, dz are expected to be in nm, no conversion needed
        for (int i = 0; i < res.atomCount; ++i) {
            model::MCAtom& atom = state.atoms[res.atomStart + i];
            atom.x += dx;  // All coordinates are in nm
            atom.y += dy;
            atom.z += dz;
            applyPBC(state, atom.x, atom.y, atom.z);
        }
        updateGeometricCenter(state, res);
    }
}

void MCCore::applyPBC(const model::MCState& state, float& x, float& y, float& z) const {
    // Box dimensions and coordinates are in nm, no conversion needed
    x -= state.info.box[0] * std::floor(x / state.info.box[0]);
    y -= state.info.box[1] * std::floor(y / state.info.box[1]);
    z -= state.info.box[2] * std::floor(z / state.info.box[2]);
}

void MCCore::updateGeometricCenter(const model::MCState& state, model::MCResidue& res) const {
    res.center[0] = res.center[1] = res.center[2] = 0.0f;
    
    // All coordinates are already in nm, no conversion needed
    for (int i = 0; i < res.atomCount; ++i) {
        const model::MCAtom& atom = state.atoms[res.atomStart + i];
        res.center[0] += atom.x;
        res.center[1] += atom.y;
        res.center[2] += atom.z;
    }
    
    if (res.atomCount > 0) {
        float invCount = 1.0f / res.atomCount;
        res.center[0] *= invCount;
        res.center[1] *= invCount;
        res.center[2] *= invCount;
    }
}

} // namespace montecarlo
} // namespace system
} // namespace pygcmc 