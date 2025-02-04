#include "MonteCarloSystem.hpp"
#include "system/system.hpp"
#include <cmath>
#include <algorithm>
#include <cctype>

namespace pygcmc {
namespace system {

using System = pygcmc::system::System;
using LogLevel = pygcmc::system::LogLevel;

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
        const auto& topRes = molecular->topology_residues[i];
        
        model::MCResidue mcRes;
        mcRes.atomStart = atomStart;
        mcRes.atomCount = molRes->atom_count();
        mcRes.active = true;
        
        // Convert atoms for this residue
        const auto& molAtoms = molRes->get_atoms();
        for (size_t j = 0; j < molAtoms.size(); j++) {
            const auto& molAtom = molAtoms[j];
            const auto& topAtom = molecular->topology_atoms[topRes.atoms[j]];
            
            model::MCAtom mcAtom;
            mcAtom.x = molAtom->get_x();
            mcAtom.y = molAtom->get_y();
            mcAtom.z = molAtom->get_z();
            mcAtom.charge = topAtom.charge;  // 使用 topology 中的电荷
            mcAtom.type = typeMaps.getOrAddType(topAtom.type);  // 使用 topology 中的类型
            
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

void MonteCarloSystem::addMovementMolecules(const std::vector<MovementMolecularInfo>& molecules)
{
    // --------------------------------------------------------------------
    // [1] 首先先收集"新分子"中的 residue types 和 atom types
    // --------------------------------------------------------------------
    pygcmc::model::TypeMaps preResidueTypes;
    pygcmc::model::TypeMaps preAtomTypes;

    // 先把所有新的分子类型加进去，这样它们的 index 会排在前面
    for (const auto& info : molecules) {
        if (!info.molecular) {
            throw std::runtime_error("Movement molecular data is null");
        }
        const auto& molRes = info.molecular->residues[0];
        std::string rawName = molRes->get_resname();
        std::string nameTrimmed = trim(rawName);
        std::string resName = nameTrimmed;
        std::transform(resName.begin(), resName.end(), resName.begin(),
                       [](unsigned char c){ return std::toupper(c); });

        // 将该新 residue 名字放到新的 map 里
        preResidueTypes.getOrAddType(resName);

        // 将该新 residue 的 atom type 放到新的 map 里
        const auto& topology_atoms = info.molecular->topology_atoms;
        for (const auto& top_atom : topology_atoms) {
            preAtomTypes.getOrAddType(top_atom.type);
        }

        System::log(LogLevel::DEBUG, "Expected movement residue name for molecule: '", 
                   resName, "'");
        System::log(LogLevel::DEBUG, "Movement molecule residues:\n  resname='", 
                   resName, "'");
    }

    // --------------------------------------------------------------------
    // [2] 再把原系统（state）里已经有的 residue types 和 atom types 放进来
    // --------------------------------------------------------------------
    for (const auto& t : state.residueTypes.atomTypes) {
        preResidueTypes.getOrAddType(t);
    }
    for (const auto& t : state.atomTypes.atomTypes) {
        preAtomTypes.getOrAddType(t);
    }

    // --------------------------------------------------------------------
    // [3] 新建一个新的 MCState，用上述新的 type maps
    // --------------------------------------------------------------------
    pygcmc::model::MCState newState;
    newState.info        = state.info;
    newState.forcefield  = state.forcefield;
    newState.residueTypes = preResidueTypes;  // 用我们收集好的"新顺序" residue type
    newState.atomTypes    = preAtomTypes;     // 用我们收集好的"新顺序" atom type

    // --------------------------------------------------------------------
    // [4] 先将旧系统的所有 residue/atom 重新索引到新 type maps
    //     并先保存到临时数据结构 oldResidues/oldResidueAtoms
    // --------------------------------------------------------------------
    std::vector<pygcmc::model::MCResidue> oldResidues;
    oldResidues.reserve(state.activeResidueCount);
    std::vector<std::vector<pygcmc::model::MCAtom>> oldResidueAtoms;
    oldResidueAtoms.reserve(state.activeResidueCount);

    // 调试输出：先打印旧系统的 residue type 及残基
    System::log(LogLevel::DEBUG, "\nOriginal (old) active Residues before reindex:");
    for (int i = 0; i < state.activeResidueCount; i++) {
        const auto& oldRes = state.residues[i];
        if (!oldRes.active) {
            System::log(LogLevel::DEBUG, "  Residue ", i, " inactive, skipping.");
            continue;
        }
        std::string oldTypeName = state.residueTypes.getTypeName(oldRes.type);
        System::log(LogLevel::DEBUG, "  Residue ", i,
                   " old type index ", oldRes.type,
                   " => type name '", oldTypeName, "'");
    }

    for (int i = 0; i < state.activeResidueCount; i++) {
        const auto& oldRes = state.residues[i];
        // 如果不想拷贝 inactive，可以自行判断。但你原先代码是一起拷过去再处理的
        // 这里就不分是否 active。
        // 取得旧名字
        std::string oldResName = state.residueTypes.getTypeName(oldRes.type);
        std::string oldResNameTrimmed = trim(oldResName);
        std::string oldResNameUpper   = oldResNameTrimmed;
        std::transform(
            oldResNameUpper.begin(), oldResNameUpper.end(), oldResNameUpper.begin(),
            [](unsigned char c){ return std::toupper(c); }
        );

        // 构造一个新的 Residue，里面 type 改成在 newState.residueTypes 下的 index
        pygcmc::model::MCResidue newRes = oldRes;
        newRes.type = newState.residueTypes.getOrAddType(oldResNameUpper);

        // 把它压到临时数组
        oldResidues.push_back(newRes);

        // 下面收集它的 atoms，并重新映射 atom type
        std::vector<pygcmc::model::MCAtom> theseAtoms;
        theseAtoms.reserve(newRes.atomCount);
        for (int j = 0; j < oldRes.atomCount; j++) {
            auto oldAtom = state.atoms[oldRes.atomStart + j];
            // 去旧系统里拿 atom type 字符串
            std::string oldAtomTypeName = typeMaps.getTypeName(oldAtom.type);
            // 在新的 typeMaps 下找对应的新索引
            int newAtomTypeIdx = newState.atomTypes.getOrAddType(oldAtomTypeName);
            oldAtom.type = newAtomTypeIdx;
            theseAtoms.push_back(oldAtom);
        }
        oldResidueAtoms.push_back(std::move(theseAtoms));
    }

    // --------------------------------------------------------------------
    // [5] 原函数的逻辑：把旧系统的活跃残基分组到 matchingResidues 里，
    //                  没匹配的丢到 otherResidues 里
    //    注意：我们这时候要用 oldResidues / oldResidueAtoms 来做"旧系统"的来源
    // --------------------------------------------------------------------
    // 取出运动分子对应的残基名称（先 trim，再转大写）
    std::vector<std::string> insertionResNames;
    insertionResNames.reserve(molecules.size());

    for (const auto& info : molecules) {
        if (!info.molecular) {
            throw std::runtime_error("Movement molecular data is null");
        }
        const auto& molRes = info.molecular->residues[0];
        std::string rawName = molRes->get_resname();
        std::string nameTrimmed = trim(rawName);
        std::string resName = nameTrimmed;
        std::transform(resName.begin(), resName.end(), resName.begin(),
                       [](unsigned char c){ return std::toupper(c); });
        insertionResNames.push_back(resName);
    }

    // 用来收集和分组
    std::vector<std::vector<pygcmc::model::MCResidue>> matchingResidues(molecules.size());
    std::vector<std::vector<std::vector<pygcmc::model::MCAtom>>> matchingAtoms(molecules.size());
    std::vector<pygcmc::model::MCResidue> otherResidues;
    std::vector<std::vector<pygcmc::model::MCAtom>> otherAtoms;

    // 第一遍：遍历 oldResidues，把 active 的 residue 做匹配
    for (int i = 0; i < static_cast<int>(oldResidues.size()); i++) {
        const auto& oldRes = oldResidues[i];
        if (!oldRes.active) {
            System::log(LogLevel::DEBUG, "Skipping inactive residue ", i);
            continue;
        }

        std::string oldResNameUpper = newState.residueTypes.getTypeName(oldRes.type);
        // 这里我们之前已经在上面转把 name 转成大写并 getOrAddType 了
        // 故此处 oldResNameUpper 应该已经是大写，但为了保持调试输出一致，还是留着
        System::log(LogLevel::DEBUG, "Processing residue ", i, 
                   ": upper='", oldResNameUpper, "'");

        bool found = false;
        for (size_t m = 0; m < molecules.size(); m++) {
            System::log(LogLevel::DEBUG, "   Comparing with insertionResNames[", m, "]: '", 
                       insertionResNames[m], "'");
            if (oldResNameUpper == insertionResNames[m]) {
                System::log(LogLevel::DEBUG, "   Residue ", i, " matched movement molecule index ", m);
                matchingResidues[m].push_back(oldRes);
                std::vector<pygcmc::model::MCAtom> atoms;
                atoms.reserve(oldRes.atomCount);
                for (int j = 0; j < oldRes.atomCount; j++) {
                    atoms.push_back(oldResidueAtoms[i][j]);
                }
                matchingAtoms[m].push_back(atoms);
                found = true;
                break;
            }
        }
        if (!found) {
            System::log(LogLevel::DEBUG, "Residue ", i, " did not match any movement molecule; adding to others.");
            otherResidues.push_back(oldRes);
            std::vector<pygcmc::model::MCAtom> atoms;
            atoms.reserve(oldRes.atomCount);
            for (int j = 0; j < oldRes.atomCount; j++) {
                atoms.push_back(oldResidueAtoms[i][j]);
            }
            otherAtoms.push_back(std::move(atoms));
        }
    }

    // 输出每种运动分子匹配到的残基数目
    for (size_t m = 0; m < molecules.size(); m++) {
        System::log(LogLevel::DEBUG, "Movement molecule index ", m, " ('", 
                   insertionResNames[m], "') collected ", 
                   matchingResidues[m].size(), " active residues.");
    }

    // --------------------------------------------------------------------
    // [6] 第二遍：构造 newState 的 residues 和 atoms
    // --------------------------------------------------------------------
    newState.atoms.reserve(state.atoms.size());
    newState.residues.reserve(state.residues.size());
    
    int newAtomStart = 0;
    int newResIdx = 0;

    // 先把未匹配到的残基（otherResidues）放进 newState
    for (size_t i = 0; i < otherResidues.size(); i++) {
        pygcmc::model::MCResidue newRes = otherResidues[i];
        newRes.atomStart = newAtomStart;
        const auto& atoms = otherAtoms[i];
        for (const auto& atom : atoms) {
            newState.atoms.push_back(atom);
        }
        newAtomStart += static_cast<int>(atoms.size());
        newState.residues.push_back(newRes);
        newResIdx++;
    }

    // 再把各个 molecules 匹配到的旧残基，以及新的不活跃拷贝加进去
    for (size_t m = 0; m < molecules.size(); m++) {
        const auto& molInfo = molecules[m];
        const auto& matches = matchingResidues[m];
        const auto& matchAtoms = matchingAtoms[m];

        int startIndexForThisGroup = newResIdx;
        int activeCountForThisGroup = static_cast<int>(matches.size());

        // 添加"旧系统"里匹配到的活跃残基
        for (size_t r = 0; r < matches.size(); r++) {
            pygcmc::model::MCResidue newRes = matches[r];
            newRes.atomStart = newAtomStart;
            newRes.active = true;
            newRes.fixed  = false;
            const auto& atoms = matchAtoms[r];
            for (const auto& atom : atoms) {
                newState.atoms.push_back(atom);
            }
            newAtomStart += static_cast<int>(atoms.size());
            newState.residues.push_back(newRes);
            newResIdx++;
        }

        // 添加不活跃拷贝
        const auto& molRes = molInfo.molecular->residues[0];
        const auto& molAtoms = molRes->get_atoms();
        const auto& topRes = molInfo.molecular->topology_residues[0];  // 获取第一个残基的 topology
        int atomsPerResidue = static_cast<int>(molAtoms.size());
        
        for (int c = 0; c < molInfo.maxCopies; c++) {
            pygcmc::model::MCResidue newRes;
            newRes.atomStart = newAtomStart;
            newRes.atomCount = atomsPerResidue;
            newRes.active = false;
            newRes.fixed  = false;
            // 新分子的 residue type 索引
            std::string rawName = molRes->get_resname();
            std::string nameTrimmed = trim(rawName);
            std::string resName = nameTrimmed;
            std::transform(
                resName.begin(), resName.end(), resName.begin(),
                [](unsigned char c){ return std::toupper(c); }
            );
            newRes.type = newState.residueTypes.getOrAddType(resName);

            // 新分子的 atom type 索引
            for (size_t i = 0; i < molAtoms.size(); i++) {
                const auto& molAtom = molAtoms[i];
                const auto& topAtom = molInfo.molecular->topology_atoms[topRes.atoms[i]];
                
                pygcmc::model::MCAtom mcAtom;
                mcAtom.x = molAtom->get_x();
                mcAtom.y = molAtom->get_y();
                mcAtom.z = molAtom->get_z();
                mcAtom.charge = topAtom.charge;  // 使用 topology 中的电荷
                mcAtom.type = newState.atomTypes.getOrAddType(topAtom.type);  // 使用 topology 中的类型
                newState.atoms.push_back(mcAtom);
            }
            newAtomStart += atomsPerResidue;
            newState.residues.push_back(newRes);
            newResIdx++;
        }

        // 记录在 movementResidues 里的信息
        pygcmc::model::MCMovementResidueInfo moveInfo;
        moveInfo.startIndex = startIndexForThisGroup;
        moveInfo.activeCount = activeCountForThisGroup;
        moveInfo.totalCount = activeCountForThisGroup + molInfo.maxCopies;
        moveInfo.resName = molRes->get_resname();
        newState.movementResidues.push_back(moveInfo);
    }

    // 更新 newState 的活跃计数
    newState.activeResidueCount = newResIdx;
    newState.activeAtomCount    = newAtomStart;

    // 检查容量
    if (newState.activeResidueCount > newState.info.maxResidues ||
        newState.activeAtomCount > newState.info.maxAtoms) {
        throw std::runtime_error("New state exceeds max capacity after movement insertion");
    }

    // --------------------------------------------------------------------
    // [7] 打印最终映射和残基信息
    // --------------------------------------------------------------------
    System::log(LogLevel::DEBUG, "\nFinal ResidueTypes mapping (newState):");
    for (size_t idx = 0; idx < newState.residueTypes.atomTypes.size(); idx++) {
        System::log(LogLevel::DEBUG, "  index=", idx, 
                   " name='", newState.residueTypes.atomTypes[idx], "'");
    }

    System::log(LogLevel::DEBUG, "\nFinal Residues (newState):");
    for (int i = 0; i < newState.activeResidueCount; i++) {
        const pygcmc::model::MCResidue& newRes = newState.residues[i];
        std::string typeName = newState.residueTypes.getTypeName(newRes.type);
        System::log(LogLevel::DEBUG, "  Residue ", i,
                   " has type index ", newRes.type,
                   " => type name '", typeName, "'",
                   (newRes.active ? " (ACTIVE)" : " (INACTIVE)"));
    }

    // 打印 movement residues 信息
    System::log(LogLevel::DEBUG, "\nMovement Residues Info (newState):");
    for (const auto& info : newState.movementResidues) {
        System::log(LogLevel::DEBUG, "  Movement group: name='", info.resName,
                   "' start=", info.startIndex,
                   " active=", info.activeCount,
                   " total=", info.totalCount);
    }

    // 最后，将 newState 替换进当前对象
    state = std::move(newState);
}


} // namespace system
} // namespace pygcmc 