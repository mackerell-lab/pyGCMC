#include "MonteCarloSystem.hpp"
#include <cmath>
#include <algorithm>
#include <cctype>

namespace pygcmc {
namespace system {

void MonteCarloSystem::addInitialResidues(const model::MCResidue* resVec, int resCount,
                                         const model::MCAtom* atomVec, int atomCount) {
    if (resCount > state.info.maxResidues || atomCount > state.info.maxAtoms) {
        throw std::runtime_error("Initial system exceeds max capacity");
    }

    std::copy(resVec, resVec + resCount, state.residues.begin());
    std::copy(atomVec, atomVec + atomCount, state.atoms.begin());

    state.activeResidueCount = resCount;
    state.activeAtomCount = atomCount;
}

int MonteCarloSystem::insertResidue(const model::MCResidue& res, const model::MCAtom* atoms) {
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

void MonteCarloSystem::translateResidue(int resIdx, float dx, float dy, float dz) {
    if (resIdx >= 0 && resIdx < state.activeResidueCount) {
        model::MCResidue& res = state.residues[resIdx];
        for (int i = 0; i < res.atomCount; ++i) {
            model::MCAtom& atom = state.atoms[res.atomStart + i];
            atom.x += dx;
            atom.y += dy;
            atom.z += dz;
            applyPBC(atom.x, atom.y, atom.z);
        }
        updateGeometricCenter(res);
    }
}

float MonteCarloSystem::calcNonBondedEnergy(const model::MCResidue& res1, const model::MCResidue& res2) const {
    (void)res1;  // Suppress unused parameter warning
    (void)res2;  // Suppress unused parameter warning
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

void MonteCarloSystem::updateGeometricCenter(model::MCResidue& res) {
    res.center[0] = res.center[1] = res.center[2] = 0.0f;
    
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

void MonteCarloSystem::initializeFromMolecular(const std::shared_ptr<model::Molecular>& molecular) {
    if (!molecular) {
        throw std::runtime_error("MolecularSystem has no molecular data");
    }
    
    // Set box dimensions from molecular system
    state.info.box[0] = molecular->boxDimensions[0];
    state.info.box[1] = molecular->boxDimensions[1];
    state.info.box[2] = molecular->boxDimensions[2];
    state.info.volume = state.info.box[0] * state.info.box[1] * state.info.box[2];

    // Convert residues and atoms
    std::vector<model::MCResidue> tempResidues;
    std::vector<model::MCAtom> tempAtoms;
    
    size_t atomStart = 0;
    const size_t numResidues = molecular->get_num_residues();
    
    for (size_t i = 0; i < numResidues; ++i) {
        const auto& molRes = molecular->residues[i];
        
        model::MCResidue mcRes;
        mcRes.atomStart = atomStart;
        mcRes.atomCount = molRes->atom_count();
        mcRes.active = true;
        
        // Convert atoms for this residue
        const auto& molAtoms = molRes->get_atoms();
        for (const auto& molAtom : molAtoms) {
            model::MCAtom mcAtom;
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

static inline std::string trim(const std::string &s) {
    auto start = s.begin();
    while (start != s.end() && std::isspace(*start)) {
        ++start;
    }
    auto end = s.end();
    while (end != start && std::isspace(*(end - 1))) {
        --end;
    }
    return std::string(start, end);
}
void MonteCarloSystem::addMovementMolecules(const std::vector<MovementMolecularInfo>& molecules) {
    // 打印初始状态的 residue types 映射
    std::cerr << "\n[DEBUG EXTRA] Initial ResidueTypes mapping:" << std::endl;
    for (size_t idx = 0; idx < state.residueTypes.atomTypes.size(); idx++) {
        std::cerr << "  index=" << idx 
                  << " name='" << state.residueTypes.atomTypes[idx] << "'" << std::endl;
    }

    // 打印初始状态的所有残基及其类型
    std::cerr << "\n[DEBUG EXTRA] Initial Residues:" << std::endl;
    for (int i = 0; i < state.activeResidueCount; i++) {
        const model::MCResidue& oldRes = state.residues[i];
        std::string typeName = state.residueTypes.getTypeName(oldRes.type);
        std::cerr << "  Residue " << i
                  << " has type index " << oldRes.type
                  << " => type name '" << typeName << "'"
                  << (oldRes.active ? " (ACTIVE)" : " (INACTIVE)") << std::endl;
    }

    // 创建一个新的状态来构造最终结果
    model::MCState newState;
    newState.info = state.info;
    newState.forcefield = state.forcefield;
    newState.atomTypes = state.atomTypes;
    newState.residueTypes = state.residueTypes;

    // 取出运动分子对应的残基名称（先 trim，再转大写）
    std::vector<std::string> insertionResNames;
    insertionResNames.reserve(molecules.size());
    
    // 预先将所有新的分子类型添加到类型映射中
    for (const auto& info : molecules) {
        if (!info.molecular) {
            throw std::runtime_error("Movement molecular data is null");
        }
        // 添加残基类型
        const auto& molRes = info.molecular->residues[0];
        std::string rawName = molRes->get_resname();
        std::string nameTrimmed = trim(rawName);
        std::string resName = nameTrimmed;
        std::transform(resName.begin(), resName.end(), resName.begin(),
                      [](unsigned char c){ return std::toupper(c); });
        insertionResNames.push_back(resName);
        
        // 确保残基类型已添加到映射中
        newState.residueTypes.getOrAddType(resName);
        
        // 添加原子类型
        const auto& molAtoms = molRes->get_atoms();
        for (const auto& molAtom : molAtoms) {
            typeMaps.getOrAddType(molAtom->get_type());
        }
        
        std::cerr << "[DEBUG] Expected movement residue name for molecule: '" << resName << "'" << std::endl;
        std::cerr << "[DEBUG EXTRA] Movement molecule residues:" << std::endl;
        std::cerr << "  resname='" << resName << "'" << std::endl;
    }

    // 第一遍：遍历基础系统中所有活跃残基，将它们按照类型分组
    std::vector<std::vector<model::MCResidue>> matchingResidues(molecules.size());
    std::vector<std::vector<std::vector<model::MCAtom>>> matchingAtoms(molecules.size());
    std::vector<model::MCResidue> otherResidues;
    std::vector<std::vector<model::MCAtom>> otherAtoms;

    for (int i = 0; i < state.activeResidueCount; i++) {
        const model::MCResidue& oldRes = state.residues[i];
        if (!oldRes.active) {
            std::cerr << "[DEBUG EXTRA] Skipping inactive residue " << i << std::endl;
            continue;
        }

        std::string rawOldResName = state.residueTypes.getTypeName(oldRes.type);
        std::string oldResNameTrimmed = trim(rawOldResName);
        std::string oldResNameUpper = oldResNameTrimmed;
        std::transform(oldResNameUpper.begin(), oldResNameUpper.end(), oldResNameUpper.begin(),
                       [](unsigned char c){ return std::toupper(c); });

        std::cerr << "[DEBUG EXTRA] Processing residue " << i 
                  << ": raw='" << rawOldResName 
                  << "', trimmed='" << oldResNameTrimmed
                  << "', upper='" << oldResNameUpper << "'" << std::endl;

        // 判断该残基是否匹配任一运动分子
        bool found = false;
        for (size_t m = 0; m < molecules.size(); m++) {
            std::cerr << "[DEBUG]   Comparing with insertionResNames[" << m << "]: '" 
                      << insertionResNames[m] << "'" << std::endl;
            if (oldResNameUpper == insertionResNames[m]) {
                std::cerr << "[DEBUG]   Residue " << i << " matched movement molecule index " << m << std::endl;
                matchingResidues[m].push_back(oldRes);
                std::vector<model::MCAtom> atoms;
                atoms.reserve(oldRes.atomCount);
                for (int j = 0; j < oldRes.atomCount; j++) {
                    atoms.push_back(state.atoms[oldRes.atomStart + j]);
                }
                matchingAtoms[m].push_back(atoms);
                found = true;
                break;
            }
        }
        if (!found) {
            std::cerr << "[DEBUG] Residue " << i << " did not match any movement molecule; adding to others." << std::endl;
            otherResidues.push_back(oldRes);
            std::vector<model::MCAtom> atoms;
            atoms.reserve(oldRes.atomCount);
            for (int j = 0; j < oldRes.atomCount; j++) {
                atoms.push_back(state.atoms[oldRes.atomStart + j]);
            }
            otherAtoms.push_back(atoms);
        }
    }

    // 输出每种运动分子匹配到的残基数目
    for (size_t m = 0; m < molecules.size(); m++) {
        std::cerr << "[DEBUG] Movement molecule index " << m << " ('" << insertionResNames[m]
                  << "') collected " << matchingResidues[m].size() << " active residues." << std::endl;
    }

    // 第二遍：构造新状态
    newState.atoms.reserve(state.atoms.size());
    newState.residues.reserve(state.residues.size());
    
    int newAtomStart = 0;
    int newResIdx = 0;

    // 先添加未匹配的（非运动）残基
    for (size_t i = 0; i < otherResidues.size(); i++) {
        model::MCResidue newRes = otherResidues[i];
        newRes.atomStart = newAtomStart;
        const auto& atoms = otherAtoms[i];
        for (const auto& atom : atoms) {
            newState.atoms.push_back(atom);
        }
        newAtomStart += static_cast<int>(atoms.size());
        newState.residues.push_back(newRes);
        newResIdx++;
    }

    // 再处理每种运动分子类型
    for (size_t m = 0; m < molecules.size(); m++) {
        const auto& molInfo = molecules[m];
        const auto& matches = matchingResidues[m];
        const auto& matchAtoms = matchingAtoms[m];

        int startIndexForThisGroup = newResIdx;
        int activeCountForThisGroup = static_cast<int>(matches.size());

        // 添加从基础系统中匹配到的活跃残基
        for (size_t r = 0; r < matches.size(); r++) {
            model::MCResidue newRes = matches[r];
            newRes.atomStart = newAtomStart;
            newRes.active = true;
            newRes.fixed = false;
            const auto& atoms = matchAtoms[r];
            for (const auto& atom : atoms) {
                newState.atoms.push_back(atom);
            }
            newAtomStart += static_cast<int>(atoms.size());
            newState.residues.push_back(newRes);
            newResIdx++;
        }

        // 添加不活跃的拷贝
        const auto& molRes = molInfo.molecular->residues[0];
        const auto& molAtoms = molRes->get_atoms();
        int atomsPerResidue = static_cast<int>(molAtoms.size());
        
        for (int c = 0; c < molInfo.maxCopies; c++) {
            model::MCResidue newRes;
            newRes.atomStart = newAtomStart;
            newRes.atomCount = atomsPerResidue;
            newRes.active = false;
            newRes.fixed = false;
            newRes.type = state.residueTypes.getOrAddType(molRes->get_resname());
            for (const auto& molAtom : molAtoms) {
                model::MCAtom mcAtom;
                mcAtom.x = molAtom->get_x();
                mcAtom.y = molAtom->get_y();
                mcAtom.z = molAtom->get_z();
                mcAtom.charge = molAtom->get_charge();
                mcAtom.type = typeMaps.getOrAddType(molAtom->get_type());
                newState.atoms.push_back(mcAtom);
            }
            newAtomStart += atomsPerResidue;
            newState.residues.push_back(newRes);
            newResIdx++;
        }

        model::MCMovementResidueInfo moveInfo;
        moveInfo.startIndex = startIndexForThisGroup;
        moveInfo.activeCount = activeCountForThisGroup;
        moveInfo.totalCount = activeCountForThisGroup + molInfo.maxCopies;
        moveInfo.resName = molRes->get_resname();
        newState.movementResidues.push_back(moveInfo);
    }

    newState.activeResidueCount = newResIdx;
    newState.activeAtomCount = newAtomStart;

    if (newState.activeResidueCount > newState.info.maxResidues ||
        newState.activeAtomCount > newState.info.maxAtoms) {
        throw std::runtime_error("New state exceeds max capacity after movement insertion");
    }

    // 在构建新状态后，打印最终的映射和残基信息
    std::cerr << "\n[DEBUG EXTRA] Final ResidueTypes mapping:" << std::endl;
    for (size_t idx = 0; idx < newState.residueTypes.atomTypes.size(); idx++) {
        std::cerr << "  index=" << idx 
                  << " name='" << newState.residueTypes.atomTypes[idx] << "'" << std::endl;
    }

    std::cerr << "\n[DEBUG EXTRA] Final Residues:" << std::endl;
    for (int i = 0; i < newState.activeResidueCount; i++) {
        const model::MCResidue& newRes = newState.residues[i];
        std::string typeName = newState.residueTypes.getTypeName(newRes.type);
        std::cerr << "  Residue " << i
                  << " has type index " << newRes.type
                  << " => type name '" << typeName << "'"
                  << (newRes.active ? " (ACTIVE)" : " (INACTIVE)") << std::endl;
    }

    // 打印 movement residues 信息
    std::cerr << "\n[DEBUG EXTRA] Movement Residues Info:" << std::endl;
    for (const auto& info : newState.movementResidues) {
        std::cerr << "  Movement group: name='" << info.resName
                  << "' start=" << info.startIndex
                  << " active=" << info.activeCount
                  << " total=" << info.totalCount << std::endl;
    }

    state = std::move(newState);
}


} // namespace system
} // namespace pygcmc 