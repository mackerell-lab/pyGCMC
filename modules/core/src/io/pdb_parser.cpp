// modules/core/src/io/pdb_parser.cpp

#include "pygcmc/core/io/pdb_parser.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include "pygcmc/core/utils.hpp"
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <optional>

namespace pygcmc {
namespace core {
namespace io {

std::pair<std::optional<std::vector<double>>, std::vector<IOResidue>> PDBParser::parse(const std::string& filename) {
    std::optional<std::vector<double>> cryst;
    std::vector<IOResidue> residues;
    IOResidue* current_residue = nullptr;

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw ParserError("无法打开文件: " + filename);
    }

    std::string line;
    while (std::getline(infile, line)) {
        // Parse crystal information
        if (line.substr(0, 6) == "CRYST1") {
            std::vector<double> box_dims;
            if (parse_cryst1_line(line, box_dims)) {
                cryst = box_dims;  // Store the box dimensions if successfully parsed
            }
            continue;
        }

        // Only parse ATOM and HETATM records
        if ((line.substr(0, 6) == "ATOM  " || line.substr(0, 6) == "HETATM") && line.length() > 54) {
            PDBAtom atom;
            if (parse_atom_line(line, atom) && atom.is_valid()) {
                // Check if a new residue needs to be created
                if (current_residue == nullptr || 
                    current_residue->name != atom.residue ||
                    current_residue->sequence_number != atom.sequence ||
                    current_residue->chain_id != atom.chain) {
                    
                    residues.emplace_back(atom.residue, atom.sequence, atom.chain);
                    current_residue = &residues.back();
                }

                // Add the PDBAtom to the current residue
                current_residue->atoms.push_back(atom);
            }
        }
    }

    infile.close();

    // Validate residue structure
    if (!validate_pdb_structure(residues)) {
        throw ParserError("PDB 文件的残基结构验证失败: " + filename);
    }

    return {cryst, residues};
}

bool PDBParser::parse_cryst1_line(const std::string& line, std::vector<double>& cell_params) {
    try {
        // CRYST1 format:
        // Columns  Data
        // 1-6      "CRYST1"
        // 7-15     a (Angstroms)
        // 16-24    b (Angstroms)
        // 25-33    c (Angstroms)
        // 34-40    alpha (degrees)
        // 41-47    beta (degrees)
        // 48-54    gamma (degrees)
        // 56-66    Space group
        // 67-70    Z value
        
        if (line.length() < 33) {  // Minimum length for a, b, c values
            return false;
        }

        double a = std::stod(utils::trim(line.substr(6, 9)));
        double b = std::stod(utils::trim(line.substr(15, 9)));
        double c = std::stod(utils::trim(line.substr(24, 9)));
        
        // Also parse angles if available
        double alpha = 90.0, beta = 90.0, gamma = 90.0;
        if (line.length() >= 54) {
            alpha = std::stod(utils::trim(line.substr(33, 7)));
            beta = std::stod(utils::trim(line.substr(40, 7)));
            gamma = std::stod(utils::trim(line.substr(47, 7)));
        }

        // Check for valid box dimensions and angles
        if (a <= 0.0 || b <= 0.0 || c <= 0.0 ||
            alpha <= 0.0 || alpha >= 180.0 ||
            beta <= 0.0 || beta >= 180.0 ||
            gamma <= 0.0 || gamma >= 180.0) {
            return false;
        }

        cell_params = {a, b, c, alpha, beta, gamma};
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

        std::string residue = utils::trim(line.substr(17, 4));
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
        // std::cout << "Parsed Atom - Serial: " << atom.serial
        //           << ", Name: " << atom.name
        //           << ", Residue: " << atom.residue
        //           << ", Sequence: " << atom.sequence
        //           << ", Chain: " << atom.chain
        //           << ", X: " << atom.x << ", Y: " << atom.y << ", Z: " << atom.z
        //           << ", Occupancy: " << atom.occupancy
        //           << ", Temp Factor: " << atom.temp_factor
        //           << ", Element: " << atom.element
        //           << ", Charge: " << atom.charge
        //           << ", Type: " << atom.type << std::endl;

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
            std::vector<std::pair<int, char>> seq_vec(sequences.begin(), sequences.end());
            
            std::sort(seq_vec.begin(), seq_vec.end(), 
                      [](const std::pair<int, char>& a, const std::pair<int, char>& b) -> bool {
                          if (a.first != b.first) return a.first < b.first;
                          return a.second < b.second;
                      });

            // Only check for duplicate sequence numbers and insertion codes
            for (size_t i = 1; i < seq_vec.size(); ++i) {
                const auto& prev = seq_vec[i-1];
                const auto& curr = seq_vec[i];
                
                if (prev.first == curr.first && prev.second == curr.second) {
                    return false;
                }
                // Removed gap check to allow for non-consecutive sequence numbers
            }
        }
    }
    
    return true;
}

bool PDBParser::validate_pdb_structure(const std::vector<IOResidue>& residues) {
    if (residues.empty()) {
        return false;  // Empty structure is invalid
    }

    // Track unique identifiers and sequences
    std::unordered_set<int> serials;
    std::unordered_map<char, std::map<std::string, std::set<std::pair<int, char>>>> chain_residues;
    // Format: chain -> residue_name -> set of (sequence, insertion_code)

    for (const auto& residue : residues) {
        for (const auto& atom : residue.atoms) {
            // Basic atom validation
            if (!validate_atom(atom)) {
                return false;
            }

            // Check for duplicate serial numbers
            if (!serials.insert(atom.serial).second) {
                return false;
            }

            // Track residue information
            auto& residue_map = chain_residues[residue.chain_id];
            residue_map[residue.name].insert({residue.sequence_number, residue.atoms[0].insertion_code});
        }
    }

    // Validate chain and residue organization
    return validate_chain_structure(chain_residues);
}

} // namespace io
} // namespace core
} // namespace pygcmc
