// modules/core/src/io/psf_parser.cpp

#include "pygcmc/core/io/parser.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <algorithm>

namespace pygcmc {
namespace core {
namespace io {

PSFTopology PSFParser::parse(const std::string& filename) {
    PSFTopology psf;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("无法打开文件: " + filename);
    }

    std::string line;
    bool in_bonds_section = false;

    while (std::getline(infile, line)) {
        // 检查是否进入BONDS部分
        if (line.find("BONDS") != std::string::npos) {
            in_bonds_section = true;
            continue;
        }

        if (in_bonds_section) {
            // PSF文件中的拓扑信息通常以数字开头
            if (std::isdigit(line[0])) {
                std::istringstream iss(line);
                int atom1, atom2;
                iss >> atom1 >> atom2;
                psf.bonds.emplace_back(PSFBond{atom1, atom2});
            } else {
                // 结束BONDS部分
                break;
            }
        }
    }

    infile.close();
    return psf;
}

} // namespace io
} // namespace core
} // namespace pygcmc
