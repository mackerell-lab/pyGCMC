// modules/core/src/io/psf_parser.cpp

#include "pygcmc/core/io/psf_parser.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>

namespace pygcmc {
namespace core {
namespace io {

PSFTopology PSFParser::parse(const std::string& filename) {
    PSFTopology psf;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw FileError("无法打开文件: " + filename);
    }

    std::string line;
    bool in_atoms_section = false;
    bool in_bonds_section = false;
    int expected_atoms = 0;
    int expected_bonds = 0;

    while (std::getline(infile, line)) {
        // 跳过空行
        if (line.empty()) continue;

        // 检查节段标记
        if (line.find("!NATOM") != std::string::npos) {
            std::istringstream iss(line);
            iss >> expected_atoms;
            in_atoms_section = true;
            in_bonds_section = false;
            continue;
        }
        else if (line.find("!NBOND") != std::string::npos) {
            std::istringstream iss(line);
            iss >> expected_bonds;
            in_atoms_section = false;
            in_bonds_section = true;
            continue;
        }

        // 解析原子信息
        if (in_atoms_section) {
            std::istringstream iss(line);
            PSFAtom atom;
            std::string segment, residue_name;  // 暂时不使用的字段

            // PSF格式：ID SEGID RESID RESNAME ATOMNAME TYPE CHARGE MASS
            if (!(iss >> atom.id >> segment >> atom.residue_id >> residue_name 
                     >> atom.name >> atom.type >> atom.charge >> atom.mass)) {
                continue;  // 跳过无法解析的行
            }

            if (atom.is_valid()) {
                psf.atoms.push_back(atom);
            }
        }
        // 解析键信息
        else if (in_bonds_section) {
            std::istringstream iss(line);
            int atom1, atom2;
            
            // PSF bonds格式：每行可能包含多个键对
            while (iss >> atom1 >> atom2) {
                PSFBond bond(atom1, atom2);
                if (bond.is_valid()) {
                    psf.bonds.push_back(bond);
                }
            }
        }
    }

    infile.close();

    // 验证解析结果
    if (expected_atoms > 0 && psf.atoms.size() != static_cast<size_t>(expected_atoms)) {
        throw FormatError("原子数量不匹配: 期望 " + std::to_string(expected_atoms) + 
                         ", 实际 " + std::to_string(psf.atoms.size()));
    }

    if (expected_bonds > 0 && psf.bonds.size() != static_cast<size_t>(expected_bonds)) {
        throw FormatError("键数量不匹配: 期望 " + std::to_string(expected_bonds) + 
                         ", 实际 " + std::to_string(psf.bonds.size()));
    }

    return psf;
}

} // namespace io
} // namespace core
} // namespace pygcmc
