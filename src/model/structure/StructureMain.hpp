// src/model/structure.hpp

#pragma once
#ifndef PYGCMC_MODEL_STRUCTURE_MAIN_HPP
#define PYGCMC_MODEL_STRUCTURE_MAIN_HPP

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "model/atom/AtomMain.hpp"
#include "model/residue/ResidueMain.hpp"

namespace pygcmc {
namespace model {
namespace structure {

/**
 * @brief Molecular structure class, used to represent the complete structural information of a molecule
 */
class Structure {
public:
    // Secondary structure information
    struct SecondaryStructure {
        std::string id;          // Structure identifier
        std::string initResName; // Initial residue name
        char initChainId;        // Initial chain ID
        int initSeqNum;         // Initial sequence number
        char initICode;         // Initial insertion code
        std::string endResName;  // Terminal residue name
        char endChainId;        // Terminal chain ID
        int endSeqNum;          // Terminal sequence number
        char endICode;          // Terminal insertion code
        int structureClass;     // Structure type
    };

    // Chain termination information
    struct TerminalInfo {
        char chainId;           // Chain ID
        int resSeq;            // Residue sequence number
        char iCode;            // Insertion code
        std::string resName;    // Residue name
    };

    Structure() = default;

    // Getters
    const std::vector<std::shared_ptr<atom::Atom>>& get_atoms() const { return atoms_; }
    const std::vector<std::shared_ptr<residue::Residue>>& get_residues() const { return residues_; }
    const std::vector<TerminalInfo>& get_terminals() const { return terminals_; }
    const std::map<std::string, std::vector<SecondaryStructure>>& get_helices() const { return helices_; }
    const std::map<std::string, std::vector<std::string>>& get_sheets() const { return sheets_; }
    const std::vector<std::string>& get_ssbonds() const { return ssbonds_; }
    const std::vector<double>& get_box_dimensions() const { return boxDimensions_; }

    // Setters
    void add_atom(const std::shared_ptr<atom::Atom>& atom) { atoms_.push_back(atom); }
    void add_residue(const std::shared_ptr<residue::Residue>& residue) { residues_.push_back(residue); }
    void add_terminal(const TerminalInfo& terminal) { terminals_.push_back(terminal); }
    void add_helix(const std::string& chainId, const SecondaryStructure& helix) {
        helices_[chainId].push_back(helix);
    }
    void add_sheet(const std::string& chainId, const std::string& sheetInfo) {
        sheets_[chainId].push_back(sheetInfo);
    }
    void add_ssbond(const std::string& ssbond) { ssbonds_.push_back(ssbond); }
    void set_box_dimensions(const std::vector<double>& dimensions) { boxDimensions_ = dimensions; }

    // Clear all data
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
    std::vector<std::shared_ptr<atom::Atom>> atoms_;
    std::vector<std::shared_ptr<residue::Residue>> residues_;
    std::vector<TerminalInfo> terminals_;
    std::map<std::string, std::vector<SecondaryStructure>> helices_;
    std::map<std::string, std::vector<std::string>> sheets_;
    std::vector<std::string> ssbonds_;
    std::vector<double> boxDimensions_;
};

} // namespace structure
} // namespace model
} // namespace pygcmc

#endif // PYGCMC_MODEL_STRUCTURE_MAIN_HPP
