// modules/core/src/io/pdb_parser.cpp

#include "pygcmc/core/io/pdb_parser.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>

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

        // 粒子类型可以根据原子名称或其他规则设定
        int type = 0; // 这里简单设为0，实际应用中可根据需求调整

        atom = PDBAtom(serial, name, residue, sequence, x, y, z, charge, type);
        return true;
    } catch (...) {
        return false;
    }
}

bool PDBParser::validate_pdb_structure(const std::vector<PDBAtom>& atoms) {
    // 添加验证逻辑，例如检查是否有重复的原子序号等
    return true;
}

} // namespace io
} // namespace core
} // namespace pygcmc
