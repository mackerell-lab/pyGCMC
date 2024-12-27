// modules/core/src/io/pdb_parser.cpp

#include "pygcmc/core/io/pdb_parser.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include "pygcmc/core/utils.hpp"
#include <unordered_set>
#include <unordered_map>

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
        int serial = std::stoi(line.substr(6, 5));
        std::string name = line.substr(12, 4);
        name.erase(std::remove_if(name.begin(), name.end(), ::isspace), name.end());

        std::string residue = line.substr(17, 3);
        residue.erase(std::remove_if(residue.begin(), residue.end(), ::isspace), residue.end());

        int sequence = std::stoi(line.substr(22, 4));

        double x = std::stod(line.substr(30, 8));
        double y = std::stod(line.substr(38, 8));
        double z = std::stod(line.substr(46, 8));

        // 可选字段：电荷
        double charge = 0.0;
        if (line.length() >= 80) {
            std::string chargeStr = line.substr(78, 2);
            chargeStr.erase(std::remove_if(chargeStr.begin(), chargeStr.end(), ::isspace), chargeStr.end());
            if (!chargeStr.empty()) {
                charge = std::stod(chargeStr);
            }
        }

        // 解析 type 作为字符串，从第 77-78 列提取
        std::string type = "";
        if (line.length() >= 78) {
            type = pygcmc::core::utils::trim(line.substr(76, 2));
        }

        atom = PDBAtom(serial, name, residue, sequence, x, y, z, charge, type);
        return true;
    } catch (...) {
        return false;
    }
}

bool PDBParser::validate_pdb_structure(const std::vector<PDBAtom>& atoms) {
    if (atoms.empty()) {
        return false;  // Empty structure is invalid
    }

    std::unordered_set<int> serials;
    std::unordered_map<std::string, std::unordered_set<int>> residue_sequences;

    for (const auto& atom : atoms) {
        // Check if atom is valid
        if (!atom.is_valid()) {
            return false;
        }

        // Check for duplicate serial numbers
        if (!serials.insert(atom.serial).second) {
            return false;
        }

        // Check coordinates are finite
        if (!std::isfinite(atom.x) || !std::isfinite(atom.y) || !std::isfinite(atom.z)) {
            return false;
        }

        // Track residue sequence numbers for each residue name
        residue_sequences[atom.residue].insert(atom.sequence);
    }

    // Check residue sequence continuity
    for (const auto& [residue, sequences] : residue_sequences) {
        std::vector<int> seq_nums(sequences.begin(), sequences.end());
        std::sort(seq_nums.begin(), seq_nums.end());
        
        // Check for gaps in sequence numbers
        for (size_t i = 1; i < seq_nums.size(); ++i) {
            if (seq_nums[i] - seq_nums[i-1] > 1) {
                return false;  // Gap detected in residue sequence
            }
        }
    }

    return true;
}

} // namespace io
} // namespace core
} // namespace pygcmc
