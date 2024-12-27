// modules/core/src/io/ff_parser.cpp

#include "pygcmc/core/io/ff_parser.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <algorithm>
#include <utility>
#include "pygcmc/core/string.hpp" // 保持包含路径

namespace pygcmc {
namespace core {
namespace io {

std::pair<NBMap, NBFixMap> FFParser::parse(const std::string& filename) {
    NBMap nb_dict;
    NBFixMap nbfix_dict;

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw FileError("无法打开文件: " + filename);
    }

    std::string line;
    bool in_atomtypes_section = false;
    bool in_nonbond_params_section = false;
    bool in_pairtypes_section = false;

    while (std::getline(infile, line)) {
        // 忽略注释
        size_t comment_pos = line.find(';');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // 去除行首尾空白
        line = pygcmc::core::utils::trim(line); // 如果 `trim` 在 `utils.hpp`

        if (line.empty()) {
            continue;
        }

        // 检查进入不同部分
        if (line.find("[ atomtypes ]") != std::string::npos) {
            in_atomtypes_section = true;
            in_nonbond_params_section = false;
            in_pairtypes_section = false;
            continue;
        }
        if (line.find("[ nonbond_params ]") != std::string::npos) {
            in_nonbond_params_section = true;
            in_atomtypes_section = false;
            in_pairtypes_section = false;
            continue;
        }
        if (line.find("[ pairtypes ]") != std::string::npos) {
            in_pairtypes_section = true;
            in_atomtypes_section = false;
            in_nonbond_params_section = false;
            continue;
        }

        if (in_atomtypes_section) {
            if (line.empty() || line[0] == '[') {
                in_atomtypes_section = false;
                continue;
            }

            if (!parse_atomtypes_line(line, nb_dict)) {
                throw FormatError("解析 atomtypes 行失败: " + line);
            }
        }

        if (in_nonbond_params_section) {
            if (line.empty() || line[0] == '[') {
                in_nonbond_params_section = false;
                continue;
            }

            if (!parse_nonbond_params_line(line, nbfix_dict)) {
                throw FormatError("解析 nonbond_params 行失败: " + line);
            }
        }

        if (in_pairtypes_section) {
            if (line.empty() || line[0] == '[') {
                in_pairtypes_section = false;
                continue;
            }

            if (!parse_pairtypes_line(line, nbfix_dict)) {
                throw FormatError("解析 pairtypes 行失败: " + line);
            }
        }
    }

    infile.close();

    return {nb_dict, nbfix_dict};
}

bool FFParser::parse_atomtypes_line(const std::string& line, NBMap& nb_dict) {
    std::istringstream iss(line);
    std::string name;
    int type;
    double charge, mass;
    double sigma, epsilon;

    iss >> name >> type >> charge >> mass >> sigma >> epsilon;

    if (name.empty() || iss.fail()) {
        return false;
    }

    nb_dict[name] = ForceFieldPair(sigma, epsilon);
    return true;
}

bool FFParser::parse_nonbond_params_line(const std::string& line, NBFixMap& nbfix_dict) {
    std::istringstream iss(line);
    std::string type1, type2;
    double sigma, epsilon;
    std::string dummy;  // For any additional fields

    iss >> type1 >> type2 >> sigma >> epsilon;

    if (type1.empty() || type2.empty() || iss.fail()) {
        return false;
    }

    nbfix_dict[{type1, type2}] = ForceFieldPair(sigma, epsilon);
    return true;
}

bool FFParser::parse_pairtypes_line(const std::string& line, NBFixMap& nbfix_dict) {
    // pairtypes 的解析与 nonbond_params 类似
    return parse_nonbond_params_line(line, nbfix_dict);
}

} // namespace io
} // namespace core
} // namespace pygcmc
