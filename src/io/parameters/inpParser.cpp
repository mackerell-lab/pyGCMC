// src/io/inpParser.cpp

#include "inpParser.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <regex>

namespace pygcmc {
namespace io {

// Helper function to split string by delimiter
std::vector<std::string> split(const std::string& str, char delim = ' ') {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(str);
    while (std::getline(tokenStream, token, delim)) {
        if (!token.empty()) {
            tokens.push_back(token);
        }
    }
    return tokens;
}

// Helper function to trim whitespace
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

// Helper function to parse array of floats
std::array<float, 3> parse_float_array(const std::string& str) {
    std::array<float, 3> result = {0.0f, 0.0f, 0.0f};
    std::istringstream iss(str);
    for (int i = 0; i < 3; ++i) {
        if (!(iss >> result[i])) {
            throw std::runtime_error("Failed to parse float array: " + str);
        }
    }
    return result;
}

// Helper function to parse vector of strings
std::vector<std::string> parse_string_vector(const std::string& str) {
    std::vector<std::string> result;
    std::istringstream iss(str);
    std::string item;
    while (iss >> item) {
        result.push_back(item);
    }
    return result;
}

// Helper function to parse vector of floats
std::vector<float> parse_float_vector(const std::string& str) {
    std::vector<float> result;
    std::istringstream iss(str);
    float value;
    while (iss >> value) {
        result.push_back(value);
    }
    return result;
}

void INPParser::parse_to_param(const std::string& filename, model::Param& param) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open input file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto tokens = split(line, ':');
        if (tokens.size() != 2) continue;

        std::string key = trim(tokens[0]);
        std::string value = trim(tokens[1]);

        try {
            parse_line(key, value, param);
        } catch (const std::exception& e) {
            throw std::runtime_error("Error parsing line '" + line + "': " + e.what());
        }
    }

    // Post-process and validate parameters
    validate_parameters(param);
}

void INPParser::parse_string_to_param(const std::string& content, model::Param& param) {
    std::istringstream iss(content);
    std::string line;
    while (std::getline(iss, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto tokens = split(line, ':');
        if (tokens.size() != 2) continue;

        std::string key = trim(tokens[0]);
        std::string value = trim(tokens[1]);

        try {
            parse_line(key, value, param);
        } catch (const std::exception& e) {
            throw std::runtime_error("Error parsing line '" + line + "': " + e.what());
        }
    }

    // Post-process and validate parameters
    validate_parameters(param);
}

void INPParser::parse_line(const std::string& key, const std::string& value, model::Param& param) {
    auto& file_info = param.get_file_info();
    auto& space_info = param.get_space_info();
    auto& fragment_info = param.get_fragment_info();
    auto& mc_info = param.get_mc_info();
    auto& bias_info = param.get_bias_info();
    auto& basic_info = param.get_basic_info();

    // File paths
    if (key == "par") {
        file_info.par_files.push_back(value);
    } else if (key == "fragitp") {
        file_info.fragment_top_files.push_back(value);
    } else if (key == "atomtypes") {
        file_info.atomtype_file = value;
    } else if (key == "monomerdir") {
        file_info.monomer_dir = value;
    } else if (key == "top") {
        file_info.topology_file = value;
    } else if (key == "pdb") {
        file_info.input_pdb_file = value;
    } else if (key == "protitp") {
        file_info.protein_top_files.push_back(value);
    } else if (key == "op_top") {
        file_info.output_top_file = value;
    } else if (key == "op_pdb") {
        file_info.output_pdb_file = value;
    }
    // Space parameters
    else if (key == "grid_dx") {
        space_info.grid_spacing = std::stof(value);
    } else if (key == "box_size") {
        space_info.box_size = parse_float_array(value);
        space_info.volume = space_info.box_size[0] * space_info.box_size[1] * space_info.box_size[2];
    } else if (key == "cutoff") {
        space_info.cutoff = std::stof(value);
    } else if (key == "gc_center") {
        space_info.gc_center = parse_float_array(value);
    } else if (key == "sys_center") {
        space_info.sys_center = parse_float_array(value);
    }
    // Fragment parameters
    else if (key == "fragname") {
        file_info.fragment_names = parse_string_vector(value);
    } else if (key == "fragconc") {
        fragment_info.conc_list = parse_float_vector(value);
    } else if (key == "fragmuex") {
        fragment_info.muex_list = parse_float_vector(value);
    }
    // MC parameters
    else if (key == "nprint") {
        mc_info.print_freq = std::stoi(value);
    } else if (key == "mcsteps") {
        mc_info.mc_steps = std::stoi(value);
    }
    // Bias parameters
    else if (key == "use_cavity_bias") {
        bias_info.use_cavity_bias = (value == "yes");
    } else if (key == "use_conf_bias") {
        bias_info.use_conf_bias = (value == "yes");
    }
    // Basic info parameters
    else if (key == "initcycle") {
        basic_info.init_cycle = (value == "yes");
    } else if (key == "conserve_frags") {
        basic_info.conserve_fragments = (value == "yes");
    } else if (key == "map_generation") {
        file_info.generate_maps = (value == "yes");
    } else if (key == "map_filename_prefix") {
        file_info.map_prefix = value;
    }
}

void INPParser::validate_parameters(model::Param& param) {
    auto& file_info = param.get_file_info();
    auto& fragment_info = param.get_fragment_info();
    auto& space_info = param.get_space_info();

    // Check required file paths
    if (file_info.topology_file.empty()) {
        throw std::runtime_error("Missing required parameter: top");
    }
    if (file_info.input_pdb_file.empty()) {
        throw std::runtime_error("Missing required parameter: pdb");
    }

    // Check fragment parameters consistency
    if (file_info.fragment_names.size() != fragment_info.conc_list.size()) {
        throw std::runtime_error("Inconsistent fragment parameters: fragname and fragconc sizes don't match");
    }
    if (file_info.fragment_names.size() != fragment_info.muex_list.size()) {
        throw std::runtime_error("Inconsistent fragment parameters: fragname and fragmuex sizes don't match");
    }

    // Check space parameters
    if (space_info.grid_spacing <= 0.0f) {
        throw std::runtime_error("Invalid grid_dx: must be positive");
    }
    if (space_info.cutoff <= 0.0f) {
        throw std::runtime_error("Invalid cutoff: must be positive");
    }
    if (space_info.volume <= 0.0f) {
        throw std::runtime_error("Invalid box_size: volume must be positive");
    }
}

model::Param INPParser::parse_file(const std::string& filename) {
    model::Param param;
    parse_to_param(filename, param);
    return param;
}

model::Param INPParser::parse_string(const std::string& content) {
    model::Param param;
    parse_string_to_param(content, param);
    return param;
}

} // namespace io
} // namespace pygcmc
