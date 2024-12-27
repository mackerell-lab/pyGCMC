// modules/core/src/io/pdb_parser.cpp

#include "pygcmc/core/io/pdb_parser.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include "pygcmc/core/utils.hpp"
#include <unordered_set>
#include <unordered_map>
#include <iostream> // 添加调试输出

namespace pygcmc {
namespace core {
namespace io {

std::pair<std::vector<double>, std::vector<PDBAtom>> PDBParser::parse(const std::string& filename) {
    std::vector<double> cryst;
    std::vector<PDBAtom> atoms;

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw ParserError("无法打开文件: " + filename);
    }

    std::string line;
    while (std::getline(infile, line)) {
        // 解析晶胞信息
        if (line.substr(0, 6) == "CRYST1") {
            if (!parse_cryst1_line(line, cryst)) {
                throw ParserError("解析晶胞信息失败: " + filename);
            }
            continue;
        }

        // 仅解析ATOM和HETATM记录
        if ((line.substr(0, 6) == "ATOM  " || line.substr(0, 6) == "HETATM") && line.length() > 54) {
            PDBAtom atom;
            if (parse_atom_line(line, atom) && atom.is_valid()) {
                atoms.emplace_back(atom);
                // 调试输出
                std::cout << "Parsed Atom: Serial " << atom.serial << ", Name " << atom.name 
                          << ", Type " << atom.type << std::endl;
            }
        }
    }

    infile.close();

    // 如果晶胞信息未找到，计算基于原子坐标的晶胞
    if (cryst.empty()) {
        if (atoms.empty()) {
            throw ParserError("没有找到晶胞信息且原子列表为空: " + filename);
        }
        double min_x = atoms[0].x, max_x = atoms[0].x;
        double min_y = atoms[0].y, max_y = atoms[0].y;
        double min_z = atoms[0].z, max_z = atoms[0].z;

        for (const auto& atom : atoms) {
            min_x = std::min(min_x, atom.x);
            max_x = std::max(max_x, atom.x);
            min_y = std::min(min_y, atom.y);
            max_y = std::max(max_y, atom.y);
            min_z = std::min(min_z, atom.z);
            max_z = std::max(max_z, atom.z);
        }

        cryst = {max_x - min_x, max_y - min_y, max_z - min_z};
    }

    return {cryst, atoms};
}

bool PDBParser::parse_cryst1_line(const std::string& line, std::vector<double>& cell_params) {
    try {
        double a = std::stod(line.substr(6, 9));
        double b = std::stod(line.substr(15, 9));
        double c = std::stod(line.substr(24, 9));
        cell_params = {a, b, c};
        return true;
    } catch (...) {
        return false;
    }
}

bool PDBParser::parse_atom_line(const std::string& line, PDBAtom& atom) {
    try {
        // Minimum length check (through column 54 for essential fields)
        if (line.length() < 54) {
            return false;
        }

        // Parse mandatory fields
        int serial = std::stoi(utils::trim(line.substr(6, 5)));

        std::string name = utils::trim(line.substr(12, 4));
        if (name.empty()) {
            return false;
        }

        char alt_loc = (line.length() > 16) ? line[16] : ' ';

        std::string residue = utils::trim(line.substr(17, 3));
        if (residue.empty()) {
            return false;
        }

        char chain = (line.length() > 21) ? line[21] : ' ';
        int sequence = std::stoi(utils::trim(line.substr(22, 4)));
        char insertion_code = (line.length() > 26) ? line[26] : ' ';

        // Parse coordinates (required fields)
        double x = std::stod(utils::trim(line.substr(30, 8)));
        double y = std::stod(utils::trim(line.substr(38, 8)));
        double z = std::stod(utils::trim(line.substr(46, 8)));

        // Parse optional fields with defaults
        double occupancy = 1.0;  // Default occupancy is 1.0
        if (line.length() >= 60) {
            std::string occupancyStr = utils::trim(line.substr(54, 6));
            if (!occupancyStr.empty()) {
                occupancy = std::stod(occupancyStr);
            }
        }

        double temp_factor = 0.0;
        if (line.length() >= 66) {
            std::string tempFactorStr = utils::trim(line.substr(60, 6));
            if (!tempFactorStr.empty()) {
                temp_factor = std::stod(tempFactorStr);
            }
        }

        // Parse element and charge (optional)
        std::string element;
        if (line.length() >= 78) {
            element = utils::trim(line.substr(76, 2));
            if (element.empty()) {
                // If element is not specified, try to derive it from atom name
                element = derive_element_from_name(name);
            }
        } else {
            element = derive_element_from_name(name);
        }

        std::string charge;
        if (line.length() >= 80) {
            charge = utils::trim(line.substr(78, 2));
        }

        // 创建并验证原子，将 type 设置为 element
        atom = PDBAtom(serial, name, residue, sequence, chain, alt_loc, insertion_code,
                      x, y, z, occupancy, temp_factor, element, charge, element);

        // 调试输出
        std::cout << "Parsed Atom - Serial: " << atom.serial
                  << ", Name: " << atom.name
                  << ", Residue: " << atom.residue
                  << ", Sequence: " << atom.sequence
                  << ", Chain: " << atom.chain
                  << ", X: " << atom.x << ", Y: " << atom.y << ", Z: " << atom.z
                  << ", Occupancy: " << atom.occupancy
                  << ", Temp Factor: " << atom.temp_factor
                  << ", Element: " << atom.element
                  << ", Charge: " << atom.charge
                  << ", Type: " << atom.type << std::endl;

        return atom.is_valid();

    } catch (const std::exception& e) {
        std::cerr << "Error parsing atom line: " << e.what() << std::endl;
        return false;
    }
}

std::string PDBParser::derive_element_from_name(const std::string& name) {
    if (name.empty()) return "";

    // If name starts with a digit, element is the rest
    if (std::isdigit(name[0])) {
        if (name.length() >= 2) {
            return utils::trim(name.substr(1, 1));  // Take second character
        } else {
            return "";
        }
    }

    // Otherwise take first one or two characters based on name
    if (name.length() >= 2 && std::isupper(name[1])) {
        return utils::trim(name.substr(0, 2));  // Two-letter element
    } else {
        return utils::trim(name.substr(0, 1));  // One-letter element
    }
}

bool PDBParser::validate_atom(const PDBAtom& atom) {
    // Check atom validity
    if (!atom.is_valid()) {
        return false;
    }

    // Check coordinates
    if (!std::isfinite(atom.x) || !std::isfinite(atom.y) || !std::isfinite(atom.z)) {
        return false;
    }

    // Check occupancy range
    if (atom.occupancy < 0.0 || atom.occupancy > 1.0) {
        return false;
    }

    // Check temperature factor
    if (!std::isfinite(atom.temp_factor) || atom.temp_factor < 0.0) {
        return false;
    }

    return true;
}

bool PDBParser::validate_chain_structure(const std::unordered_map<char, 
    std::map<std::string, std::set<std::pair<int, char>>>>& chain_residues) {
        
    for (const auto& [chain, residues] : chain_residues) {
        for (const auto& [residue_name, sequences] : residues) {
            // Convert sequence/insertion pairs to vector for analysis
            std::vector<std::pair<int, char>> seq_vec(sequences.begin(), sequences.end());
            
            // Sort based on sequence number and insertion code
            std::sort(seq_vec.begin(), seq_vec.end(), 
                      [](const std::pair<int, char>& a, const std::pair<int, char>& b) -> bool {
                          if (a.first != b.first) return a.first < b.first;
                          return a.second < b.second;
                      });

            // Check for sequence continuity
            for (size_t i = 1; i < seq_vec.size(); ++i) {
                const auto& prev = seq_vec[i-1];
                const auto& curr = seq_vec[i];
                
                // If same sequence number, must have different insertion codes
                if (prev.first == curr.first && prev.second == curr.second) {
                    return false;
                }
                
                // Check for unreasonable gaps (more than 1 residue)
                if (curr.first - prev.first > 1) {
                    return false;  // Gap detected in residue sequence
                }

                // Optional: Handle insertion codes appropriately if needed
            }
        }
    }
    
    return true;
}

bool PDBParser::validate_pdb_structure(const std::vector<PDBAtom>& atoms) {
    if (atoms.empty()) {
        return false;  // Empty structure is invalid
    }

    // Track unique identifiers and sequences
    std::unordered_set<int> serials;
    std::unordered_map<char, std::map<std::string, std::set<std::pair<int, char>>>> chain_residues;
    // Format: chain -> residue_name -> set of (sequence, insertion_code)

    for (const auto& atom : atoms) {
        // Basic atom validation
        if (!validate_atom(atom)) {
            return false;
        }

        // Check for duplicate serial numbers
        if (!serials.insert(atom.serial).second) {
            return false;
        }

        // Track residue information
        auto& residue_map = chain_residues[atom.chain];
        residue_map[atom.residue].insert({atom.sequence, atom.insertion_code});
    }

    // Validate chain and residue organization
    return validate_chain_structure(chain_residues);
}

} // namespace io
} // namespace core
} // namespace pygcmc
