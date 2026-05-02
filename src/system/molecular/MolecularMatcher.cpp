#include "MolecularMatcher.hpp"
#include "MolecularValidator.hpp"
#include "../log/LogMain.hpp"
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <queue>

namespace pygcmc {
namespace system {
namespace molecular {

using LogMain = pygcmc::system::log::LogMain;
using LogLevel = pygcmc::system::common::LogLevel;

// Static member definition
const std::set<std::string> MolecularMatcher::standard_amino_acids_ = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "HSE", "HSP", "HSC",  // Histidine different protonation states
    "CYX",  // Disulfide bond form cysteine
    "HID", "HIE", "HIP"   // Histidine in CHARMM force field
};

MolecularMatcher::MolecularMatcher()
    : validator_(std::make_unique<MolecularValidator>()) {
}

bool MolecularMatcher::matchResidueSequence(
    const std::vector<std::shared_ptr<model::Residue>>& pdb_residues,
    const std::shared_ptr<model::Topology>& topology,
    size_t start_idx,
    size_t& matched_count) {

    matched_count = 0;
    const size_t top_num_residues = topology->get_num_residues();

    // If the remaining residue count is insufficient, return false directly
    if (start_idx + top_num_residues > pdb_residues.size()) {
        return false;
    }

    // Check if the residue sequence matches
    for (size_t i = 0; i < top_num_residues; ++i) {
        const auto& pdb_res = pdb_residues[start_idx + i];
        const auto& top_res = topology->get_residue(static_cast<int>(i));

        if (pdb_res->get_resname() != top_res.name) {
            return false;
        }

        try {
            validator_->verifyAtomTypes(pdb_res, top_res, topology);
        } catch (const std::runtime_error&) {
            return false;
        }

        matched_count++;
    }

    return true;
}

void MolecularMatcher::performDetailedMatching(
    const std::shared_ptr<model::Molecular>& molecular,
    const std::shared_ptr<model::Topology>& topology) {

    const auto num_residues = static_cast<size_t>(topology->get_num_residues());

    // Build residue connection relationships
    auto residue_connections = buildResidueConnections(topology);

    // Find connected chains
    auto chains = findConnectedChains(residue_connections, num_residues);

    // Calculate amino acid count on each chain
    auto chain_aa_count = calculateChainAACount(chains, topology);

    // Establish mapping for non-protein residues
    std::map<std::string, std::vector<size_t>> mol_type_indices;
    for (size_t i = 0; i < num_residues; ++i) {
        const auto& top_res = topology->get_residue(static_cast<int>(i));
        mol_type_indices[top_res.name].push_back(i);
    }

    // Maintain current index for each type of non-protein molecule
    std::map<std::string, size_t> current_mol_index;

    // Verify atom count and type for each residue
    for (const auto& res : molecular->residues) {
        bool found_matching_res = false;

        // Find current residue's chain
        int current_res_id = -1;
        std::set<int>* current_chain = nullptr;

        for (size_t i = 0; i < num_residues; ++i) {
            const auto& top_res = topology->get_residue(static_cast<int>(i));
            if (top_res.name == res->get_resname()) {
                current_res_id = static_cast<int>(i);
                break;
            }
        }

        if (current_res_id >= 0) {
            for (auto& chain : chains) {
                if (chain.find(current_res_id) != chain.end()) {
                    current_chain = &chain;
                    break;
                }
            }
        }

        // Check if it's a protein residue
        bool is_amino_acid = isStandardAminoAcid(res->get_resname());
        bool is_protein = false;
        if (current_chain && chain_aa_count[*current_chain] > 3) {
            is_protein = is_amino_acid;
        }

        if (is_protein) {
            // For protein residues, match by name and number exactly
            for (size_t i = 0; i < num_residues; ++i) {
                const auto& top_res = topology->get_residue(static_cast<int>(i));
                if (top_res.name == res->get_resname() && top_res.number == res->get_ires()) {
                    found_matching_res = true;
                    validator_->verifyAtomTypes(res, top_res, topology);
                    break;
                }
            }
        } else {
            // For non-protein residues, match by molecule type in sequence
            const std::string& resname = res->get_resname();
            auto it = mol_type_indices.find(resname);
            if (it != mol_type_indices.end()) {
                LogMain::log(LogLevel::DEBUG, "Found ", it->second.size(), " instances of ",
                          resname, " in topology, current index: ", current_mol_index[resname]);
            } else {
                LogMain::log(LogLevel::DEBUG, "No instances of ", resname, " found in topology");
            }

            if (it != mol_type_indices.end() && current_mol_index[resname] < it->second.size()) {
                size_t top_res_idx = it->second[current_mol_index[resname]];
                const auto& top_res = topology->get_residue(static_cast<int>(top_res_idx));
                found_matching_res = true;
                current_mol_index[resname]++;

                validator_->verifyAtomTypes(res, top_res, topology);
            }
        }

        if (!found_matching_res) {
            std::stringstream ss;
            ss << "Could not find matching residue in topology for "
               << res->get_resname() << " " << res->get_ires();
            throw std::runtime_error(ss.str());
        }
    }
}

std::map<int, std::set<int>> MolecularMatcher::buildResidueConnections(
    const std::shared_ptr<model::Topology>& topology) {

    std::map<int, std::set<int>> residue_connections;
    const auto& bonds = topology->get_bonds();

    for (const auto& bond : bonds) {
        const auto& atom1 = topology->get_atom(bond.atom1);
        const auto& atom2 = topology->get_atom(bond.atom2);
        if (atom1.residue_id != atom2.residue_id) {
            residue_connections[atom1.residue_id].insert(atom2.residue_id);
            residue_connections[atom2.residue_id].insert(atom1.residue_id);
        }
    }

    return residue_connections;
}

std::vector<std::set<int>> MolecularMatcher::findConnectedChains(
    const std::map<int, std::set<int>>& residue_connections,
    size_t num_residues) {

    std::set<int> visited_residues;
    std::vector<std::set<int>> chains;

    for (size_t i = 0; i < num_residues; ++i) {
        int res_id = static_cast<int>(i);
        if (visited_residues.find(res_id) != visited_residues.end()) {
            continue;
        }

        // Find all residues connected to the current residue
        std::set<int> current_chain;
        std::queue<int> to_visit;
        to_visit.push(res_id);

        while (!to_visit.empty()) {
            int current = to_visit.front();
            to_visit.pop();

            if (visited_residues.find(current) != visited_residues.end()) {
                continue;
            }

            visited_residues.insert(current);
            current_chain.insert(current);

            // Add connected residues to the queue
            if (residue_connections.find(current) != residue_connections.end()) {
                for (int connected : residue_connections.at(current)) {
                    if (visited_residues.find(connected) == visited_residues.end()) {
                        to_visit.push(connected);
                    }
                }
            }
        }

        if (!current_chain.empty()) {
            chains.push_back(current_chain);
        }
    }

    return chains;
}

std::map<std::set<int>, int> MolecularMatcher::calculateChainAACount(
    const std::vector<std::set<int>>& chains,
    const std::shared_ptr<model::Topology>& topology) {

    std::map<std::set<int>, int> chain_aa_count;

    for (const auto& chain : chains) {
        int aa_count = 0;
        for (int res_id : chain) {
            const auto& res = topology->get_residue(res_id);
            if (isStandardAminoAcid(res.name)) {
                aa_count++;
            }
        }
        chain_aa_count[chain] = aa_count;
    }

    return chain_aa_count;
}

bool MolecularMatcher::isStandardAminoAcid(const std::string& resname) const {
    return standard_amino_acids_.find(resname) != standard_amino_acids_.end();
}

} // namespace molecular
} // namespace system
} // namespace pygcmc
