// modules/core/src/io/psf_parser.cpp

#include "pygcmc/core/io/psf_parser.hpp"
#include "pygcmc/core/utils.hpp" // 添加这一行以包含 utils.hpp
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
    bool in_bonds_section = false;

    while (std::getline(infile, line)) {
        // 解析并去除注释
        size_t comment_pos = line.find(';');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // 去除行首尾空白
        line = pygcmc::core::utils::trim(line); // 使用完整的命名空间前缀

        if (line.empty()) {
            continue;
        }

        // 检查是否进入 BONDS 部分
        if (line.find("BONDS") != std::string::npos) {
            in_bonds_section = true;
            continue;
        }

        if (in_bonds_section) {
            // PSF 文件中的拓扑信息通常以数字开头
            if (std::isdigit(line[0])) {
                std::istringstream iss(line);
                int atom1, atom2;
                iss >> atom1 >> atom2;
                if (atom1 > 0 && atom2 > 0 && atom1 != atom2) {
                    psf.bonds.emplace_back(PSFBond{atom1, atom2});
                }
            } else {
                // 结束 BONDS 部分
                in_bonds_section = false;
            }
        }
    }

    infile.close();

    // 验证拓扑结构
    if (!psf.is_valid()) {
        throw FormatError("解析的 PSF 拓扑结构无效");
    }

    return psf;
}

} // namespace io
} // namespace core
} // namespace pygcmc
