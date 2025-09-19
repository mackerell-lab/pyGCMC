// src/io/parameters/InpParserMain.cpp

#include "InpParserMain.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace pygcmc {
namespace io {
namespace parameters {

model::Param InpParserMain::parse_file(const std::string& filename) {
    model::Param param;
    parse_to_param(filename, param);
    return param;
}

model::Param InpParserMain::parse_string(const std::string& content) {
    model::Param param;
    parse_string_to_param(content, param);
    return param;
}

void InpParserMain::parse_to_param(const std::string& filename, model::Param& param) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open input file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        line = InpParserStructures::trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto tokens = InpParserStructures::split(line, ':');
        if (tokens.size() != 2) continue;

        std::string key = InpParserStructures::trim(tokens[0]);
        std::string value = InpParserStructures::trim(tokens[1]);

        try {
            parse_line(key, value, param);
        } catch (const std::exception& e) {
            throw std::runtime_error("Error parsing line '" + line + "': " + e.what());
        }
    }
    validate_parameters(param);
}

void InpParserMain::parse_string_to_param(const std::string& content, model::Param& param) {
    std::istringstream iss(content);
    std::string line;
    while (std::getline(iss, line)) {
        line = InpParserStructures::trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto tokens = InpParserStructures::split(line, ':');
        if (tokens.size() != 2) continue;

        std::string key = InpParserStructures::trim(tokens[0]);
        std::string value = InpParserStructures::trim(tokens[1]);

        try {
            parse_line(key, value, param);
        } catch (const std::exception& e) {
            throw std::runtime_error("Error parsing line '" + line + "': " + e.what());
        }
    }
    validate_parameters(param);
}

void InpParserMain::parse_line(const std::string& key, const std::string& value, model::Param& param) {
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
    } else if (key == "box_size" || key == "box") {
        space_info.box_size = InpParserStructures::parse_float_array(value);
        space_info.volume = space_info.box_size[0] * space_info.box_size[1] * space_info.box_size[2];
    } else if (key == "cutoff") {
        space_info.cutoff = std::stof(value);
    } else if (key == "gc_center") {
        space_info.gc_center = InpParserStructures::parse_float_array(value);
    } else if (key == "sys_center") {
        space_info.sys_center = InpParserStructures::parse_float_array(value);
    }
    // Fragment parameters
    else if (key == "fragname") {
        file_info.fragment_names = InpParserStructures::parse_string_vector(value);
    } else if (key == "fragconc") {
        fragment_info.conc_list = InpParserStructures::parse_float_vector(value);
    } else if (key == "fragmuex") {
        fragment_info.muex_list = InpParserStructures::parse_float_vector(value);
    }
    // MC parameters
    else if (key == "nprint") {
        mc_info.print_freq = std::stoi(value);
    } else if (key == "mcsteps") {
        mc_info.mc_steps = std::stoi(value);
    } else if (key == "temperature") {
        mc_info.temperature = std::stof(value);
        // Calculate beta from temperature (beta = 1/(kB*T) in kJ/mol units)
        mc_info.beta = 1.0f / (0.00831446f * mc_info.temperature);
    } else if (key == "eqsteps") {
        // Store equilibration steps if needed
        // Currently not used in MCParams, but parsed for compatibility
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

void InpParserMain::validate_parameters(model::Param& param) {
    auto& file_info = param.get_file_info();
    auto& fragment_info = param.get_fragment_info();
    auto& space_info = param.get_space_info();

    // For GCMC simulations with fragment top files (fragitp), we don't require PDB or TOP files
    // Only validate them if they're provided
    bool has_fragitp = !file_info.fragment_top_files.empty();
    
    // If no fragment files, then we need both PDB and TOP for a valid simulation
    // Check in the order expected by tests: top first, then pdb
    if (!has_fragitp) {
        if (file_info.topology_file.empty()) {
            throw std::runtime_error("Missing required parameter: top");
        }
        if (file_info.input_pdb_file.empty()) {
            throw std::runtime_error("Missing required parameter: pdb");
        }
    }
    
    // Fragment parameter consistency checks
    if (file_info.fragment_names.size() != fragment_info.conc_list.size()) {
        throw std::runtime_error("Inconsistent fragment parameters: fragname and fragconc sizes don't match");
    }
    if (file_info.fragment_names.size() != fragment_info.muex_list.size()) {
        throw std::runtime_error("Inconsistent fragment parameters: fragname and fragmuex sizes don't match");
    }
    
    // Grid and space checks
    if (space_info.grid_spacing <= 0.0f) {
        throw std::runtime_error("Invalid grid_dx: must be positive");
    }
    if (space_info.cutoff <= 0.0f) {
        throw std::runtime_error("Invalid cutoff: must be positive");
    }
    
    // Only check volume if box_size was actually provided (not all zeros)
    bool box_provided = (space_info.box_size[0] > 0.0f || 
                         space_info.box_size[1] > 0.0f || 
                         space_info.box_size[2] > 0.0f);
    if (box_provided && space_info.volume <= 0.0f) {
        throw std::runtime_error("Invalid box_size: volume must be positive");
    }
}

} // namespace parameters
} // namespace io
} // namespace pygcmc