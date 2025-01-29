// src/model/structure.hpp

#ifndef PYGCMC_MODEL_STRUCTURE_HPP
#define PYGCMC_MODEL_STRUCTURE_HPP

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "model/atom.hpp"
#include "model/residue.hpp"

namespace pygcmc {
namespace model {

/**
 * @brief 分子结构类,用于表示分子的完整结构信息
 */
class Structure {
public:
    // 二级结构信息
    struct SecondaryStructure {
        std::string id;          // 结构标识符
        std::string initResName; // 起始残基名
        char initChainId;        // 起始链ID
        int initSeqNum;         // 起始序号
        char initICode;         // 起始插入码
        std::string endResName;  // 终止残基名
        char endChainId;        // 终止链ID
        int endSeqNum;          // 终止序号
        char endICode;          // 终止插入码
        int structureClass;     // 结构类型
    };

    // 链终止信息
    struct TerminalInfo {
        char chainId;           // 链ID
        int resSeq;            // 残基序号
        char iCode;            // 插入码
        std::string resName;    // 残基名
    };

    Structure() = default;

    // Getters
    const std::vector<std::shared_ptr<Atom>>& getAtoms() const { return atoms_; }
    const std::vector<std::shared_ptr<Residue>>& getResidues() const { return residues_; }
    const std::vector<TerminalInfo>& getTerminals() const { return terminals_; }
    const std::map<std::string, std::vector<SecondaryStructure>>& getHelices() const { return helices_; }
    const std::map<std::string, std::vector<std::string>>& getSheets() const { return sheets_; }
    const std::vector<std::string>& getSSBonds() const { return ssbonds_; }
    const std::vector<double>& getBoxDimensions() const { return boxDimensions_; }

    // Setters
    void addAtom(const std::shared_ptr<Atom>& atom) { atoms_.push_back(atom); }
    void addResidue(const std::shared_ptr<Residue>& residue) { residues_.push_back(residue); }
    void addTerminal(const TerminalInfo& terminal) { terminals_.push_back(terminal); }
    void addHelix(const std::string& chainId, const SecondaryStructure& helix) { 
        helices_[chainId].push_back(helix); 
    }
    void addSheet(const std::string& chainId, const std::string& sheetInfo) {
        sheets_[chainId].push_back(sheetInfo);
    }
    void addSSBond(const std::string& ssbond) { ssbonds_.push_back(ssbond); }
    void setBoxDimensions(const std::vector<double>& dimensions) { boxDimensions_ = dimensions; }

    // 清除所有数据
    void clear() {
        atoms_.clear();
        residues_.clear();
        terminals_.clear();
        helices_.clear();
        sheets_.clear();
        ssbonds_.clear();
        boxDimensions_.clear();
    }

private:
    std::vector<std::shared_ptr<Atom>> atoms_;
    std::vector<std::shared_ptr<Residue>> residues_;
    std::vector<TerminalInfo> terminals_;
    std::map<std::string, std::vector<SecondaryStructure>> helices_;
    std::map<std::string, std::vector<std::string>> sheets_;
    std::vector<std::string> ssbonds_;
    std::vector<double> boxDimensions_;
};

} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_STRUCTURE_HPP
