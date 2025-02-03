#include "MonteCarloSystem.hpp"
#include <cmath>

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

void MonteCarloSystem::addMovementMolecules(const std::vector<MovementMolecularInfo>& molecules) {
    // 对每种移动分子类型进行处理
    for (const auto& molInfo : molecules) {
        const auto& molecular = molInfo.molecular;
        if (!molecular) {
            throw std::runtime_error("Movement molecular data is null");
        }
        
        // 1. 获取当前分子类型的残基名称
        std::string resName = molecular->residues[0]->get_resname();
        
        // 2. 用于临时保存匹配残基及其原子数据
        std::vector<model::MCResidue> movedResidues;
        std::vector<std::vector<model::MCAtom>> movedResidueAtoms;  // 每个残基对应一组原子
        
        // 3. 扫描当前活性残基（索引 [0, activeResidueCount)），
        //    对于匹配 resName 的残基，采用 swap‐and‐pop 方式删除，并保存其数据
        int i = 0;
        while (i < state.activeResidueCount) {
            if (state.residueTypes.getTypeName(state.residues[i].type) == resName) {
                // 保存当前匹配残基及对应的原子数据
                model::MCResidue removedRes = state.residues[i];
                std::vector<model::MCAtom> savedAtoms;
                savedAtoms.reserve(removedRes.atomCount);
                for (int j = 0; j < removedRes.atomCount; ++j) {
                    savedAtoms.push_back(state.atoms[removedRes.atomStart + j]);
                }
                movedResidues.push_back(removedRes);
                movedResidueAtoms.push_back(savedAtoms);
                
                // 采用 swap-and-pop 删除残基 i
                int atomStart = state.residues[i].atomStart;
                int atomCount = state.residues[i].atomCount;
                // 将活性原子区域末尾的这块原子拷贝填入当前位置（如果不在目标位置的话）
                if (atomStart != state.activeAtomCount - atomCount) {
                    for (int j = 0; j < atomCount; ++j) {
                        state.atoms[atomStart + j] = state.atoms[state.activeAtomCount - atomCount + j];
                    }
                }
                state.activeAtomCount -= atomCount;
                
                // 如果当前残基不是最后一个，则用最后一个残基来填补空缺
                if (i != state.activeResidueCount - 1) {
                    state.residues[i] = state.residues[state.activeResidueCount - 1];
                    // 注意：交换过来的残基原来的 atomStart 已经指向正确的原子块
                }
                state.activeResidueCount--;
                // 删除后当前 i 处新来的残基还未检测，因此不自增 i
            } else {
                ++i;
            }
        }
        
        // 4. 将之前保存的匹配残基追加到活性区域末尾
        //    此处，新残基将依次占用从 state.activeAtomCount 开始的原子区域
        int newResStart = state.activeResidueCount;
        int newAtomStart = state.activeAtomCount;
        for (size_t idx = 0; idx < movedResidues.size(); ++idx) {
            model::MCResidue movedRes = movedResidues[idx];
            // 重设原子起始位置
            movedRes.atomStart = newAtomStart;
            // 为保险起见，确保该残基原子数与模板一致
            movedRes.atomCount = static_cast<int>(molecular->residues[0]->get_atoms().size());
            movedRes.fixed = false; // 标记为可移动
            movedRes.active = true; // 保持活性
            
            // 将保存的原子数据复制到新区域
            const std::vector<model::MCAtom>& savedAtoms = movedResidueAtoms[idx];
            for (int j = 0; j < movedRes.atomCount; ++j) {
                state.atoms[newAtomStart + j] = savedAtoms[j];
            }
            newAtomStart += movedRes.atomCount;
            state.residues[newResStart++] = movedRes;
        }
        
        // 5. 添加新的非活动副本（用于后续 gcmc 插入）
        int inactiveStart = newResStart;  // 非活动副本开始的索引
        int numInactiveToAdd = molInfo.maxCopies;  // 要添加的非活动副本数量
        const auto& molRes = molecular->residues[0];
        const auto& molAtoms = molRes->get_atoms();
        int atomsPerResidue = static_cast<int>(molAtoms.size());
        for (int i = 0; i < numInactiveToAdd; ++i) {
            model::MCResidue newRes;
            newRes.atomStart = newAtomStart;
            newRes.atomCount = atomsPerResidue;
            newRes.active = false;  // 非活动标记
            newRes.fixed = false;   // 非固定状态（后续可移动）
            newRes.type = state.residueTypes.getOrAddType(resName);
            
            // 复制每个原子
            for (const auto& molAtom : molAtoms) {
                model::MCAtom mcAtom;
                mcAtom.x = molAtom->get_x();
                mcAtom.y = molAtom->get_y();
                mcAtom.z = molAtom->get_z();
                mcAtom.charge = molAtom->get_charge();
                mcAtom.type = typeMaps.getOrAddType(molAtom->get_type());
                state.atoms[newAtomStart++] = mcAtom;
            }
            
            state.residues[newResStart++] = newRes;
        }
        
        // 6. 记录这种类型移动分子的相关信息
        model::MCMovementResidueInfo moveInfo;
        // movedResidues.size() 为原来匹配的活性残基数
        moveInfo.startIndex = inactiveStart - movedResidues.size();  // 活性残基的起始位置
        moveInfo.activeCount = static_cast<int>(movedResidues.size());
        moveInfo.totalCount = moveInfo.activeCount + numInactiveToAdd;  // 总计：原活性 + 新 inactive 副本
        moveInfo.resName = resName;
        state.movementResidues.push_back(moveInfo);
        
        // 7. 更新全局计数器，保证后续操作中 active 残基与原子区域连续
        state.activeAtomCount = newAtomStart;
        state.activeResidueCount = newResStart;
    }
}


} // namespace system
} // namespace pygcmc 