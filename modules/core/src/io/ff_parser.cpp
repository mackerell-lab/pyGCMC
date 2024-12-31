// modules/core/src/io/ff_parser.cpp

#include "pygcmc/core/io/ff_parser.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype>
#include <algorithm>
#include <cmath>

namespace pygcmc {
namespace core {
namespace io {

bool FFParser::parse(const std::string& filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        std::cerr << "Failed to open force field file: " << filename << std::endl;
        return false;
    }

    std::string line;
    bool inNonbondedSection = false;
    bool inNBFIXSection = false;

    while (std::getline(ifs, line)) {
        // Skip empty lines
        if (line.empty()) continue;

        // Remove comments and trim whitespace
        auto commentPos = line.find('!');
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
        }
        auto startPos = line.find_first_not_of(" \t\r\n");
        if (startPos == std::string::npos) continue;
        line.erase(0, startPos);
        auto endPos = line.find_last_not_of(" \t\r\n");
        if (endPos != std::string::npos) {
            line.erase(endPos + 1);
        }

        // Skip empty lines after removing comments
        if (line.empty()) continue;

        // Convert to uppercase for keyword comparison
        std::string uline = line;
        std::transform(uline.begin(), uline.end(), uline.begin(), ::toupper);

        // Detect sections
        if (uline.rfind("NONBONDED", 0) == 0) {
            inNonbondedSection = true;
            inNBFIXSection = false;
            continue;
        }
        else if (uline.rfind("NBFIX", 0) == 0) {
            inNonbondedSection = false;
            inNBFIXSection = true;
            continue;
        }
        else if (uline == "END") {
            inNonbondedSection = false;
            inNBFIXSection = false;
            continue;
        }

        // Parse lines based on current section
        if (inNonbondedSection) {
            parse_nonbonded_line(line);
        }
        else if (inNBFIXSection) {
            parse_nbfix_line(line);
        }
    }

    ifs.close();
    return true;
}

void FFParser::parse_nonbonded_line(const std::string& line) {
    // 跳过注释行等
    if (line.empty()) return;
    if (line[0] == '-' || line[0] == '*' || line[0] == '@' || line[0] == '#' || line[0] == '!')
        return;

    // 我们把整行拆成 token，再一个个识别
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    {
        std::string tk;
        while (iss >> tk) {
            tokens.push_back(tk);
        }
    }
    if (tokens.empty()) return;

    // 如果第一项就是 cutnb / ctofnb / ... 之类的关键字，
    // 那说明这是一个"全局非键参数"行。
    // 也可能它们全在一行：
    //   cutnb 14.0 ctofnb 12.0 ctonnb 10.0 ...
    // 我们就循环遍历 tokens 并解析。
    bool recognized_global_params = false;
    for (size_t i = 0; i < tokens.size(); i++) {
        std::string tk = tokens[i];
        // 转成小写做比较
        std::string lower;
        lower.resize(tk.size());
        std::transform(tk.begin(), tk.end(), lower.begin(),
                       [](unsigned char c){return std::tolower(c);});

        if (lower == "cutnb") {
            if (i+1 < tokens.size()) {
                cutnb_ = std::stod(tokens[++i]); // 读下一个作为数值
                recognized_global_params = true;
            }
        }
        else if (lower == "ctofnb") {
            if (i+1 < tokens.size()) {
                ctofnb_ = std::stod(tokens[++i]);
                recognized_global_params = true;
            }
        }
        else if (lower == "ctonnb") {
            if (i+1 < tokens.size()) {
                ctonnb_ = std::stod(tokens[++i]);
                recognized_global_params = true;
            }
        }
        else if (lower == "eps") {
            if (i+1 < tokens.size()) {
                eps_ = std::stod(tokens[++i]);
                recognized_global_params = true;
            }
        }
        else if (lower == "e14fac") {
            if (i+1 < tokens.size()) {
                e14fac_ = std::stod(tokens[++i]);
                recognized_global_params = true;
            }
        }
        else if (lower == "wmin") {
            if (i+1 < tokens.size()) {
                wmin_ = std::stod(tokens[++i]);
                recognized_global_params = true;
            }
        }
    }

    // 如果我们发现这行**只**是cutnb/ctofnb...这些关键字，并且都被成功解析，
    // 那就认为这是"全局非键参数"行，不必再解析原子类型；直接 return
    if (recognized_global_params && tokens.size() < 5) {
        return;
    }

    // 如果 tokens.size() >= 4，尝试当作 "atomType ignore epsilon rminHalf"
    if (tokens.size() >= 4 && !recognized_global_params) {
        // Special handling for SILCS parameters (LP and LQ)
        if (tokens[0] == "LP" || tokens[0] == "LQ") {
            // For LP and LQ, we use the NBFIX parameters
            // The nonbonded parameters are just placeholders
            nonbonded_params_[tokens[0]] = ForceFieldPair(0.0, 0.0);
            return;
        }

        std::string atomType = tokens[0];
        // 忽略第二个参数 (ignore)
        double epsilon = std::stod(tokens[2]);
        double rminHalf = std::stod(tokens[3]);

        double rmin = rminHalf * 2.0;
        nonbonded_params_[atomType] = ForceFieldPair(rmin, std::fabs(epsilon));

        // 如果有1-4参数，也可以存储
        if (tokens.size() >= 7) {
            // 目前1-4参数未使用，暂时跳过
            // tokens[4] = ignore2
            // tokens[5] = eps14
            // tokens[6] = rmin14Half
        }
    }
}

void FFParser::parse_nbfix_line(const std::string& line) {
    // Skip header lines
    if (line.find("Emin") != std::string::npos || 
        line.find("kcal/mol") != std::string::npos ||
        line.find("NBFIX") != std::string::npos) {
        return;
    }

    // Skip lines that start with special characters
    if (line[0] == '-' || line[0] == '*' || line[0] == '@' || line[0] == '#' || line[0] == '!') {
        return;
    }

    // Parse NBFIX line
    std::istringstream iss(line);
    std::string atom_type1, atom_type2;
    double epsilon, rmin;

    if (iss >> atom_type1 >> atom_type2 >> epsilon >> rmin) {
        // In CHARMM format, epsilon is negative
        auto key = std::make_pair(atom_type1, atom_type2);
        nbfix_params_[key] = ForceFieldPair(rmin, std::abs(epsilon));
        
        // Add reverse pair
        auto key_rev = std::make_pair(atom_type2, atom_type1);
        nbfix_params_[key_rev] = ForceFieldPair(rmin, std::abs(epsilon));
        
        // Debug output
        // std::cout << "Parsed NBFIX: " << atom_type1 << " - " << atom_type2 
        //           << " | epsilon=" << std::abs(epsilon) 
        //           << ", rmin=" << rmin << std::endl;
    }
}

int FFParser::update_pdb_atoms(std::vector<PDBAtom>& atoms) const {
    int updated = 0;
    for (auto& atom : atoms) {
        if (atom.topo_type.empty()) {
            continue;
        }
        auto it = nonbonded_params_.find(atom.topo_type);
        if (it != nonbonded_params_.end()) {
            // 仅在参数未设置（为 nan）时更新
            if (std::isnan(atom.forcefield_epsilon) || std::isnan(atom.forcefield_rmin)) {
                atom.forcefield_epsilon = it->second.epsilon;
                atom.forcefield_rmin = it->second.rmin;
                updated++;
            }
        }
    }
    return updated;
}

int FFParser::update_pdb_atoms(std::vector<PDBAtom*>& atoms) const {
    int updated = 0;
    for (auto* atom : atoms) {
        if (atom->topo_type.empty()) {
            continue;
        }
        auto it = nonbonded_params_.find(atom->topo_type);
        if (it != nonbonded_params_.end()) {
            // 仅在参数未设置（为 nan）时更新
            if (std::isnan(atom->forcefield_epsilon) || std::isnan(atom->forcefield_rmin)) {
                atom->forcefield_epsilon = it->second.epsilon;
                atom->forcefield_rmin = it->second.rmin;
                updated++;
            }
        }
    }
    return updated;
}

} // namespace io
} // namespace core
} // namespace pygcmc
