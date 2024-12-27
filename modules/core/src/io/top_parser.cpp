// modules/core/src/io/top_parser.cpp

#include "pygcmc/core/io/top_parser.hpp"
#include "pygcmc/core/utils.hpp" // 确保包含 utils.hpp
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <algorithm>

namespace pygcmc {
namespace core {
namespace io {

Topology TopParser::parse(const std::string& filename) {
    Topology topology;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw FileError("无法打开文件: " + filename);
    }

    std::string line;
    bool in_atomtypes_section = false;

    while (std::getline(infile, line)) {
        // 忽略注释
        size_t comment_pos = line.find(';');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // 去除行首尾空白
        line = pygcmc::core::utils::trim(line); // 使用完整命名空间前缀

        if (line.empty()) {
            continue;
        }

        // 检查是否进入 [ atomtypes ] 部分
        if (line.find("[ atomtypes ]") != std::string::npos) {
            in_atomtypes_section = true;
            continue;
        }

        if (in_atomtypes_section) {
            if (line.empty() || line[0] == '[') {
                // 结束 [ atomtypes ] 部分
                in_atomtypes_section = false;
                continue;
            }

            if (!parse_atomtypes_section(line, topology)) {
                throw FormatError("解析 atomtypes 行失败: " + line);
            }
        }
    }

    infile.close();

    // 验证拓扑结构
    if (!topology.is_valid()) {
        throw FormatError("解析的 TOP 拓扑结构无效");
    }

    return topology;
}

bool TopParser::parse_atomtypes_section(const std::string& line, Topology& top) {
    std::istringstream iss(line);
    std::string name;
    int type;
    double charge, mass;

    iss >> name >> type >> charge >> mass;

    if (name.empty()) {
        return false;
    }

    try {
        TopAtomType atom_type(name, type, charge, mass);
        if (!atom_type.is_valid()) {
            return false;
        }
        top.atom_types.emplace_back(atom_type);
        return true;
    } catch (const FormatError& e) {
    // 记录或处理无效的力场参数
    return false;
}
}

} // namespace io
} // namespace core
} // namespace pygcmc
