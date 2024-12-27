// modules/core/src/io/itp_parser.cpp

#include "pygcmc/core/io/itp_parser.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <algorithm>

namespace pygcmc {
namespace core {
namespace io {

std::vector<ITPAtom> ITPParser::parse(const std::string& filename) {
    std::vector<ITPAtom> itp_atoms;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw FileError("无法打开文件: " + filename);
    }

    std::string line;
    bool in_atoms_section = false;

    while (std::getline(infile, line)) {
        // 忽略注释
        size_t comment_pos = line.find(';');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // 去除行首尾空白
        line = utils::trim(line);

        if (line.empty()) {
            continue;
        }

        // 检查是否进入 [ atoms ] 部分
        if (line.find("[ atoms ]") != std::string::npos) {
            in_atoms_section = true;
            continue;
        }

        if (in_atoms_section) {
            if (line.empty() || line[0] == '[') {
                // 结束 [ atoms ] 部分
                in_atoms_section = false;
                continue;
            }

            std::istringstream iss(line);
            int serial, resid;
            std::string name, type, resname;
            double charge;

            iss >> serial >> name >> resid >> resname >> type >> charge;

            if (!name.empty()) {
                ITPAtom atom(name, type, resid, resname, charge);
                if (atom.is_valid()) {
                    itp_atoms.emplace_back(atom);
                }
            }
        }
    }

    infile.close();
    return itp_atoms;
}

} // namespace io
} // namespace core
} // namespace pygcmc
