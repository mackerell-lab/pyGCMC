#include "MolecularValidator.hpp"
#include "../log/LogMain.hpp"
#include <stdexcept>
#include <sstream>
#include <iostream>
#include <map>
#include <set>
#include <cctype>

namespace pygcmc {
namespace system {
namespace molecular {

using LogMain = pygcmc::system::log::LogMain;
using LogLevel = pygcmc::system::common::LogLevel;

void MolecularValidator::verifyAtomTypes(
    const std::shared_ptr<model::Residue>& pdb_res,
    const model::TopologyResidue& top_res,
    const std::shared_ptr<model::Topology>& topology) {
    
    const auto& pdb_atoms = pdb_res->get_atoms();
    if (pdb_atoms.size() != top_res.atoms.size()) {
        std::stringstream ss;
        ss << "Inconsistent number of atoms in residue " << pdb_res->get_resname()
           << " " << pdb_res->get_ires() << ": Structure has "
           << pdb_atoms.size() << " atoms, but Topology has "
           << top_res.atoms.size() << " atoms";
        throw std::runtime_error(ss.str());
    }
    
    for (size_t i = 0; i < pdb_atoms.size(); ++i) {
        const auto& pdb_atom = pdb_atoms[i];
        const auto& top_atom = topology->get_atom(static_cast<int>(top_res.atoms[i]));
        
        // Extract element from PDB atom name
        std::string pdb_element = pdb_atom->get_element();
        if (pdb_element.empty()) {
            pdb_element = pdb_atom->get_type();
        }
        
        // Extract element from topology atom type
        std::string top_element = top_atom.type;
        
        // Compare first letter (converted to uppercase)
        char pdb_first = std::toupper(pdb_element[0]);
        char top_first = std::toupper(top_element[0]);
        
        if (pdb_first != top_first) {
            std::stringstream ss;
            ss << "Mismatched atom elements in residue " << pdb_res->get_resname()
               << " " << pdb_res->get_ires() << ": Structure has "
               << pdb_first << " (from " << pdb_atom->get_type()
               << "), but Topology has " << top_first
               << " (from " << top_atom.type << ")";
            throw std::runtime_error(ss.str());
        }
    }
}

void MolecularValidator::validateCombination(
    const std::shared_ptr<model::Structure>& structure,
    const std::shared_ptr<model::Topology>& topology) {
    
    if (!structure || !topology) {
        throw std::invalid_argument("Structure and Topology pointers cannot be null");
    }

    const auto& atoms = structure->get_atoms();
    const auto& residues = structure->get_residues();
    const auto num_atoms = static_cast<size_t>(topology->get_num_atoms());
    const auto num_residues = static_cast<size_t>(topology->get_num_residues());

    // Verify total atom count
    if (atoms.size() != num_atoms) {
        throw std::runtime_error(generateMismatchError(structure, topology, "atoms"));
    }

    // Verify total residue count
    if (residues.size() != num_residues) {
        throw std::runtime_error(generateMismatchError(structure, topology, "residues"));
    }

    // Verify atom type first letter match
    for (size_t i = 0; i < num_atoms; ++i) {
        const auto& pdb_atom = atoms[i];
        const auto& top_atom = topology->get_atom(static_cast<int>(i));
        
        // Extract element from PDB atom name
        std::string pdb_element = pdb_atom->get_element();
        if (pdb_element.empty()) {
            pdb_element = pdb_atom->get_type();
        }
        
        // Extract element from topology atom type
        std::string top_element = top_atom.type;
        
        // Compare first letter (converted to uppercase)
        char pdb_first = std::toupper(pdb_element[0]);
        char top_first = std::toupper(top_element[0]);
        
        if (pdb_first != top_first) {
            std::stringstream ss;
            ss << "Mismatched atom elements at index " << i << ": Structure has " 
               << pdb_first << " (from " << pdb_atom->get_type() 
               << "), but Topology has " << top_first 
               << " (from " << top_atom.type << ")";
            throw std::runtime_error(ss.str());
        }
    }
}

void MolecularValidator::validateMultipleCombination(
    const std::shared_ptr<model::Structure>& structure,
    const std::vector<std::shared_ptr<model::Topology>>& topologies,
    size_t total_atoms,
    size_t total_residues) {
    
    if (!structure || topologies.empty()) {
        throw std::invalid_argument("Structure and Topologies cannot be null/empty");
    }

    const auto& atoms = structure->get_atoms();
    const auto& residues = structure->get_residues();

    // Verify if totals match
    if (atoms.size() != total_atoms) {
        std::stringstream ss;
        ss << "Total number of atoms mismatch: Structure has "
           << atoms.size() << " atoms, but combined topologies have "
           << total_atoms << " atoms";
        throw std::runtime_error(ss.str());
    }
    
    if (residues.size() != total_residues) {
        std::stringstream ss;
        ss << "Total number of residues mismatch: Structure has "
           << residues.size() << " residues, but combined topologies have "
           << total_residues << " residues";
        throw std::runtime_error(ss.str());
    }
}

std::string MolecularValidator::extractElement(const std::string& atom_name) const {
    if (atom_name.empty()) {
        return "";
    }
    return std::string(1, std::toupper(atom_name[0]));
}

std::string MolecularValidator::generateMismatchError(
    const std::shared_ptr<model::Structure>& structure,
    const std::shared_ptr<model::Topology>& topology,
    const std::string& error_type) const {
    
    std::stringstream ss;
    const auto& residues = structure->get_residues();
    const auto num_residues = static_cast<size_t>(topology->get_num_residues());
    
    if (error_type == "atoms") {
        const auto& atoms = structure->get_atoms();
        const auto num_atoms = static_cast<size_t>(topology->get_num_atoms());
        
        ss << "Inconsistent total number of atoms: Structure has " 
           << atoms.size() << " atoms, but Topology has " 
           << num_atoms << " atoms\n"
           << "This mismatch suggests that the structure file contains additional molecules "
           << "that are not present in the topology file.\n";
    } else if (error_type == "residues") {
        ss << "Inconsistent total number of residues: Structure has " 
           << residues.size() << " residues, but Topology has " 
           << num_residues << " residues\n"
           << "This mismatch suggests that the structure file contains additional residues "
           << "that are not present in the topology file.\n";
    }
    
    ss << "Structure residues (in PDB order):";
    
    // List residues from structure file in PDB order
    for (const auto& res : residues) {
        ss << "\n  " << res->get_resname() << " " << res->get_ires();
    }
    
    ss << "\n\nTopology residues:";
    // List residues from topology file in order
    for (size_t i = 0; i < num_residues; ++i) {
        const auto& res = topology->get_residue(static_cast<int>(i));
        ss << "\n  " << res.name << " " << res.number;
    }
    
    // Add residue statistics
    ss << "\n\nResidue count summary:";
    ss << "\nStructure:";
    std::map<std::string, int> struct_res_count;
    for (const auto& res : residues) {
        struct_res_count[res->get_resname()]++;
    }
    for (const auto& [resname, count] : struct_res_count) {
        ss << "\n  " << resname << ": " << count;
    }
    
    ss << "\nTopology:";
    std::map<std::string, int> top_res_count;
    for (size_t i = 0; i < num_residues; ++i) {
        const auto& res = topology->get_residue(static_cast<int>(i));
        top_res_count[res.name]++;
    }
    for (const auto& [resname, count] : top_res_count) {
        ss << "\n  " << resname << ": " << count;
    }
    
    return ss.str();
}

} // namespace molecular
} // namespace system
} // namespace pygcmc 