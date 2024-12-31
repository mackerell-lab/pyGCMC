// modules/core/src/io/psf_parser.cpp

#include "pygcmc/core/io/psf_parser.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>

namespace pygcmc {
namespace core {
namespace io {

bool PSFParser::parse(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return false;
    }

    std::vector<std::string> atom_lines;
    bool in_atoms_section = false;
    std::string line;
    int natoms = 0;

    while (std::getline(file, line)) {
        // Skip empty lines and comments starting with *
        if (line.empty() || line[0] == '*') {
            continue;
        }

        // Check for NATOM section
        if (line.find("!NATOM") != std::string::npos) {
            // Extract number of atoms
            std::istringstream iss(line);
            iss >> natoms;
            std::cout << "Found " << natoms << " atoms" << std::endl;
            in_atoms_section = true;
            continue;
        }

        // End of atoms section
        if (in_atoms_section && line.find("!") != std::string::npos) {
            in_atoms_section = false;
            continue;
        }

        if (in_atoms_section && !line.empty()) {
            atom_lines.push_back(line);
        }
    }

    return parse_atoms_section(atom_lines);
}

bool PSFParser::parse_atoms_section(const std::vector<std::string>& lines) {
    atoms_.clear();
    atom_index_.clear();
    id_index_.clear();

    for (const auto& line : lines) {
        std::istringstream iss(line);
        PSFAtom atom;

        // PSF format: ID SEGID RESID RESNAME ATOMNAME ATOMTYPE CHARGE MASS
        if (!(iss >> atom.id 
                  >> atom.segment 
                  >> atom.residue_number 
                  >> atom.residue 
                  >> atom.name 
                  >> atom.type 
                  >> atom.charge 
                  >> atom.mass)) 
        {
            std::cerr << "Failed to parse atom line: " << line << std::endl;
            continue;
        }

        std::cout << "Parsed atom: " << atom.residue << " " << atom.name 
                  << " type=" << atom.type 
                  << " charge=" << atom.charge 
                  << " mass=" << atom.mass << std::endl;

        // 存入容器
        size_t idx = atoms_.size();
        atoms_.push_back(atom);

        // 建立索引: residue_name -> (residue_number -> (atom_name -> idx))
        atom_index_[atom.residue][atom.residue_number][atom.name] = idx;
        
        // 建立 atom_id 索引
        id_index_[atom.id] = idx;
    }

    return !atoms_.empty();
}

bool PSFParser::get_atom_properties(const std::string& residue_name,
                                  const std::string& atom_name,
                                  double& charge,
                                  double& mass) const {
    // Find residue
    auto res_it = atom_index_.find(residue_name);
    if (res_it == atom_index_.end()) {
        std::cerr << "[get_atom_properties] Residue not found: " 
                  << residue_name << std::endl;
        return false;
    }

    // Get the first defined residue number for this residue type
    auto first_res_num_it = res_it->second.begin();
    if (first_res_num_it == res_it->second.end()) {
        std::cerr << "[get_atom_properties] No residue numbers found for: "
                  << residue_name << std::endl;
        return false;
    }

    // Find atom name
    auto atom_it = first_res_num_it->second.find(atom_name);
    if (atom_it == first_res_num_it->second.end()) {
        std::cerr << "[get_atom_properties] Atom not found: "
                  << residue_name << " " << atom_name << std::endl;
        return false;
    }

    const PSFAtom& atom = atoms_[atom_it->second];
    charge = atom.charge;
    mass = atom.mass;
    std::cout << "Found atom: " << residue_name << " " << atom_name 
              << " charge=" << charge << " mass=" << mass << std::endl;
    return true;
}

bool PSFParser::get_atom_properties(const std::string& residue_name,
                                  int residue_number,
                                  const std::string& atom_name,
                                  double& charge,
                                  double& mass) const {
    // Find residue
    auto res_it = atom_index_.find(residue_name);
    if (res_it == atom_index_.end()) {
        std::cerr << "[get_atom_properties] Residue not found: " 
                  << residue_name << std::endl;
        return false;
    }

    // Find residue number
    auto res_num_it = res_it->second.find(residue_number);
    if (res_num_it == res_it->second.end()) {
        std::cerr << "[get_atom_properties] Residue number not found: "
                  << residue_name << " " << residue_number << std::endl;
        return false;
    }

    // Find atom name
    auto atom_it = res_num_it->second.find(atom_name);
    if (atom_it == res_num_it->second.end()) {
        std::cerr << "[get_atom_properties] Atom not found: "
                  << residue_name << " " << residue_number << " " << atom_name << std::endl;
        return false;
    }

    const PSFAtom& atom = atoms_[atom_it->second];
    charge = atom.charge;
    mass = atom.mass;
    std::cout << "Found atom: " << residue_name << " " << residue_number << " " << atom_name 
              << " charge=" << charge << " mass=" << mass << std::endl;
    return true;
}

int PSFParser::update_pdb_atoms(std::vector<PDBAtom>& pdb_atoms) const {
    int updated = 0;
    for (auto& pdb_atom : pdb_atoms) {
        // 1) 先尝试通过 serial 号匹配
        if (pdb_atom.serial > 0) {
            auto it = id_index_.find(pdb_atom.serial);
            if (it != id_index_.end()) {
                const PSFAtom& psf_atom = atoms_[it->second];
                // 验证残基名和原子名是否匹配
                if (psf_atom.residue == pdb_atom.residue && 
                    psf_atom.name == pdb_atom.name) 
                {
                    pdb_atom.topo_type = psf_atom.type;
                    pdb_atom.topo_charge = psf_atom.charge;
                    pdb_atom.topo_mass = psf_atom.mass;
                    pdb_atom.chain = psf_atom.segment.empty() ? ' ' : psf_atom.segment[0];
                    updated++;
                    continue;
                }
            }
        }

        // 2) 如果通过 serial 号匹配失败，尝试通过 residue + sequence + name 匹配
        auto res_it = atom_index_.find(pdb_atom.residue);
        if (res_it == atom_index_.end()) {
            std::cerr << "Residue not found: " << pdb_atom.residue << std::endl;
            continue;
        }

        // 查找 residue_number
        auto res_num_it = res_it->second.find(pdb_atom.sequence);
        if (res_num_it == res_it->second.end()) {
            // 如果精确匹配不到 residue_number，尝试找到第一个包含该原子名的残基
            bool found = false;
            for (const auto& [res_num, atoms_map] : res_it->second) {
                auto atom_it = atoms_map.find(pdb_atom.name);
                if (atom_it != atoms_map.end()) {
                    res_num_it = res_it->second.find(res_num);
                    found = true;
                    break;
                }
            }
            if (!found) {
                std::cerr << "Residue number not found: " 
                          << pdb_atom.residue << " " << pdb_atom.sequence << std::endl;
                continue;
            }
        }

        // 查找 atom_name
        auto atom_it = res_num_it->second.find(pdb_atom.name);
        if (atom_it == res_num_it->second.end()) {
            std::cerr << "Atom not found: " << pdb_atom.residue << " " 
                      << pdb_atom.name << std::endl;
            continue;
        }

        // 取出 PSFAtom 并更新 PDBAtom
        const PSFAtom& psf_atom = atoms_[atom_it->second];
        pdb_atom.topo_type = psf_atom.type;
        pdb_atom.topo_charge = psf_atom.charge;
        pdb_atom.topo_mass = psf_atom.mass;
        pdb_atom.chain = psf_atom.segment.empty() ? ' ' : psf_atom.segment[0];

        std::cout << "Updated atom: " << pdb_atom.residue << " " << pdb_atom.name
                  << " type=" << pdb_atom.topo_type 
                  << " charge=" << pdb_atom.topo_charge
                  << " mass=" << pdb_atom.topo_mass 
                  << " chain=" << pdb_atom.chain << std::endl;
        updated++;
    }
    return updated;
}

std::map<std::string, std::set<std::string>> PSFParser::get_missing_topology_info(
    const std::vector<PDBAtom>& atoms) const {
    std::map<std::string, std::set<std::string>> missing_info;
    
    // First check which residues don't exist in PSF file
    for (const auto& atom : atoms) {
        if (atom_index_.find(atom.residue) == atom_index_.end()) {
            missing_info[atom.residue].insert(atom.name);
            continue;
        }
    }
    
    // Create a copy of the atoms to update with topology info
    std::vector<PDBAtom> atoms_copy = atoms;
    update_pdb_atoms(atoms_copy);
    
    // Check which atoms are still missing topology info
    for (size_t i = 0; i < atoms.size(); ++i) {
        if (!atoms_copy[i].has_topology_info() && 
            atom_index_.find(atoms[i].residue) != atom_index_.end()) {
            missing_info[atoms[i].residue].insert(atoms[i].name);
        }
    }
    
    return missing_info;
}

} // namespace io
} // namespace core
} // namespace pygcmc
