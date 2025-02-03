#include "MonteCarloSystem.hpp"
#include <cmath>

namespace gcmc {

void MonteCarloSystem::addInitialResidues(const Residue* resVec, int resCount,
                                         const Atom* atomVec, int atomCount) {
    if (resCount > state.info.maxResidues || atomCount > state.info.maxAtoms) {
        throw std::runtime_error("Initial system exceeds max capacity");
    }

    std::copy(resVec, resVec + resCount, state.residues.begin());
    std::copy(atomVec, atomVec + atomCount, state.atoms.begin());

    state.activeResidueCount = resCount;
    state.activeAtomCount = atomCount;
}

int MonteCarloSystem::insertResidue(const Residue& res, const Atom* atoms) {
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

bool MonteCarloSystem::removeResidue(int resIdx) {
    if (resIdx < 0 || resIdx >= state.activeResidueCount || !state.residues[resIdx].active) {
        return false;
    }

    // Get residue info and mark it as inactive
    Residue& res = state.residues[resIdx];
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

void MonteCarloSystem::translateResidue(int resIdx, float dx, float dy, float dz) {
    if (resIdx >= 0 && resIdx < state.activeResidueCount) {
        Residue& res = state.residues[resIdx];
        for (int i = 0; i < res.atomCount; ++i) {
            Atom& atom = state.atoms[res.atomStart + i];
            atom.x += dx;
            atom.y += dy;
            atom.z += dz;
            applyPBC(atom.x, atom.y, atom.z);
        }
        updateGeometricCenter(res);
    }
}

float MonteCarloSystem::calcNonBondedEnergy([[maybe_unused]] const Residue& res1, [[maybe_unused]] const Residue& res2) const {
    float energy = 0.0f;
    // TODO: Implement LJ + Coulomb with periodic boundary conditions
    return energy;
}

float MonteCarloSystem::calcTotalEnergy() const {
    float energy = 0.0f;
    // TODO: Implement total energy calculation
    return energy;
}

void MonteCarloSystem::applyPBC(float& x, float& y, float& z) const {
    x -= state.info.box[0] * std::floor(x / state.info.box[0]);
    y -= state.info.box[1] * std::floor(y / state.info.box[1]);
    z -= state.info.box[2] * std::floor(z / state.info.box[2]);
}

float MonteCarloSystem::getMinImageDistSqr(float dx, float dy, float dz) const {
    dx -= state.info.box[0] * std::round(dx / state.info.box[0]);
    dy -= state.info.box[1] * std::round(dy / state.info.box[1]);
    dz -= state.info.box[2] * std::round(dz / state.info.box[2]);
    return dx*dx + dy*dy + dz*dz;
}

void MonteCarloSystem::updateGeometricCenter(Residue& res) {
    res.center[0] = res.center[1] = res.center[2] = 0.0f;
    
    for (int i = 0; i < res.atomCount; ++i) {
        const Atom& atom = state.atoms[res.atomStart + i];
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

void MonteCarloSystem::initializeFromMolecular(const std::shared_ptr<pygcmc::model::Molecular>& molecular) {
    if (!molecular) {
        throw std::runtime_error("MolecularSystem has no molecular data");
    }
    
    // Set box dimensions from molecular system
    state.info.box[0] = molecular->boxDimensions[0];
    state.info.box[1] = molecular->boxDimensions[1];
    state.info.box[2] = molecular->boxDimensions[2];
    state.info.volume = state.info.box[0] * state.info.box[1] * state.info.box[2];

    // Convert residues and atoms
    std::vector<Residue> tempResidues;
    std::vector<Atom> tempAtoms;
    
    size_t atomStart = 0;
    const size_t numResidues = molecular->get_num_residues();
    
    for (size_t i = 0; i < numResidues; ++i) {
        const auto& molRes = molecular->residues[i];
        
        Residue mcRes;
        mcRes.atomStart = atomStart;
        mcRes.atomCount = molRes->atom_count();
        mcRes.active = true;
        
        // Convert atoms for this residue
        const auto& molAtoms = molRes->get_atoms();
        for (const auto& molAtom : molAtoms) {
            Atom mcAtom;
            mcAtom.x = molAtom->get_x();
            mcAtom.y = molAtom->get_y();
            mcAtom.z = molAtom->get_z();
            mcAtom.charge = molAtom->get_charge();
            mcAtom.type = typeMaps.getOrAddType(molAtom->get_type());
            
            tempAtoms.push_back(mcAtom);
        }

        // Set residue type
        mcRes.type = state.residueTypes.getOrAddType(molRes->get_resname());

        // Calculate center of mass
        mcRes.center[0] = mcRes.center[1] = mcRes.center[2] = 0.0f;
        for (const auto& atom : molAtoms) {
            mcRes.center[0] += atom->get_x();
            mcRes.center[1] += atom->get_y();
            mcRes.center[2] += atom->get_z();
        }
        
        if (mcRes.atomCount > 0) {
            float invCount = 1.0f / mcRes.atomCount;
            mcRes.center[0] *= invCount;
            mcRes.center[1] *= invCount;
            mcRes.center[2] *= invCount;
        }
        
        tempResidues.push_back(mcRes);
        atomStart += mcRes.atomCount;
    }

    // Check capacity
    if (tempResidues.size() > static_cast<size_t>(state.info.maxResidues) ||
        tempAtoms.size() > static_cast<size_t>(state.info.maxAtoms)) {
        throw std::runtime_error("Initial system exceeds max capacity");
    }

    // Initialize the system with converted data
    addInitialResidues(tempResidues.data(), tempResidues.size(),
                      tempAtoms.data(), tempAtoms.size());
}

} // namespace gcmc 