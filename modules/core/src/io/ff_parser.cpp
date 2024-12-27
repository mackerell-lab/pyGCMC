// modules/core/src/io/ff_parser.cpp

#include "pygcmc/core/io/parser.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <algorithm>
#include <utility>

namespace pygcmc {
namespace core {
namespace io {

std::pair<FFParser::NBMap, FFParser::NBFixMap> FFParser::parse(const std::string& filename) {
    FFParser::NBMap nb_dict;
    FFParser::NBFixMap nbfix_dict;

    std::ifstream infile(filename);
    if (!infile.is_open()) {
        throw std::runtime_error("无法打开文件: " + filename);
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

            std::istringstream iss(line);
            std::string name;
            int type;
            double charge, mass;
            double sigma, epsilon;

            iss >> name >> type >> charge >> mass >> sigma >> epsilon;

            if (!name.empty()) {
                nb_dict[name] = ForceFieldPair{sigma, epsilon};
            }
        }

        if (in_nonbond_params_section) {
            if (line.empty() || line[0] == '[') {
                in_nonbond_params_section = false;
                continue;
            }

            std::istringstream iss(line);
            std::string type1, type2;
            double sigma, epsilon, rmin;

            iss >> type1 >> type2 >> sigma >> epsilon >> rmin;

            if (!type1.empty() && !type2.empty()) {
                nbfix_dict[{type1, type2}] = ForceFieldPair{sigma, epsilon};
            }
        }

        if (in_pairtypes_section) {
            if (line.empty() || line[0] == '[') {
                in_pairtypes_section = false;
                continue;
            }

            std::istringstream iss(line);
            std::string type1, type2;
            double sigma, epsilon, rmin;

            iss >> type1 >> type2 >> sigma >> epsilon >> rmin;

            if (!type1.empty() && !type2.empty()) {
                nbfix_dict[{type1, type2}] = ForceFieldPair{sigma, epsilon};
            }
        }
    }

    infile.close();
    return {nb_dict, nbfix_dict};
}

} // namespace io
} // namespace core
} // namespace pygcmc
