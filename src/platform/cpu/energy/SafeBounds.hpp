#ifndef SAFE_BOUNDS_HPP
#define SAFE_BOUNDS_HPP

#include <algorithm>
#include <cassert>

// 安全的边界检查宏，用于快速修复越界问题

// 在Release模式下仅做边界限制，Debug模式下会assert
#ifdef NDEBUG
  #define SAFE_ATOM_INDEX(i, atoms) std::min(i, static_cast<int>(atoms.size()) - 1)
  #define SAFE_RESIDUE_INDEX(r, residues) std::min(r, static_cast<int>(residues.size()) - 1)
#else
  #define SAFE_ATOM_INDEX(i, atoms) \
    (assert(i >= 0 && i < static_cast<int>(atoms.size())), i)
  #define SAFE_RESIDUE_INDEX(r, residues) \
    (assert(r >= 0 && r < static_cast<int>(residues.size())), r)
#endif

// 安全的循环上界
#define SAFE_ATOM_LIMIT(state) \
    std::min(state.activeAtomCount, static_cast<int>(state.atoms.size()))

#define SAFE_RESIDUE_LIMIT(state) \
    std::min(state.activeResidueCount, static_cast<int>(state.residues.size()))

// 验证residue的atom范围
inline bool isResidueAtomsValid(const MCResidue& res, int atomsSize) {
    return res.atomStart >= 0 && 
           res.atomStart + res.atomCount <= atomsSize;
}

// 安全的atom循环范围
inline std::pair<int, int> safeAtomRange(const MCResidue& res, int atomsSize) {
    int start = std::max(0, res.atomStart);
    int end = std::min(res.atomStart + res.atomCount, atomsSize);
    return {start, std::max(start, end)};
}

#endif // SAFE_BOUNDS_HPP