// modules/core/src/io/itp_parser.cpp

#include "pygcmc/core/io/parser.hpp"
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
        throw std::runtime_error("无法打开文件: " + filename);
    }

    std::string line;
    bool in_atoms_section = false;

    while (std::getline(infile, line)) {
        // 忽略注释
        size_t comment_pos = line.find(';');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // 检查是否进入ATOMS部分
        if (line.find("[ atoms ]") != std::string::npos) {
            in_atoms_section = true;
            continue;
        }

        if (in_atoms_section) {
            if (line.empty() || line[0] == '[') {
                // 结束ATOMS部分
                break;
            }

            std::istringstream iss(line);
            int serial, resid;
            std::string name, type, resname;
            double charge;

            iss >> serial >> name >> resid >> resname >> type >> charge;

            if (!name.empty()) {
                itp_atoms.emplace_back(ITPAtom{name, type, resid, resname, charge});
            }
        }
    }

    infile.close();
    return itp_atoms;
}

} // namespace io
} // namespace core
} // namespace pygcmc
