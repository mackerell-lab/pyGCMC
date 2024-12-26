// modules/core/src/io/top_parser.cpp

#include "pygcmc/core/io/parser.hpp"
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
        throw std::runtime_error("无法打开文件: " + filename);
    }

    std::string line;
    bool in_atomtypes_section = false;

    while (std::getline(infile, line)) {
        // 忽略注释
        size_t comment_pos = line.find(';');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // 检查是否进入ATOMTYPES部分
        if (line.find("[ atomtypes ]") != std::string::npos) {
            in_atomtypes_section = true;
            continue;
        }

        if (in_atomtypes_section) {
            if (line.empty() || line[0] == '[') {
                // 结束ATOMTYPES部分
                break;
            }

            std::istringstream iss(line);
            std::string name;
            int type;
            double charge, mass, sigma, epsilon;

            iss >> name >> type >> charge >> mass >> sigma >> epsilon;

            if (!name.empty()) {
                topology.atom_types.emplace_back(TopAtomType{name, type, charge, mass});
            }
        }
    }

    infile.close();
    return topology;
}

} // namespace io
} // namespace core
} // namespace pygcmc
