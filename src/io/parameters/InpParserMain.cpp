// src/io/parameters/InpParserMain.cpp

#include "InpParserMain.hpp"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

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

    auto& basic_info = param.get_basic_info();
    basic_info.inp_keys_seen.clear();
    basic_info.inp_keys_handled.clear();
    basic_info.inp_keys_unknown.clear();
    basic_info.inp_keys_ignored.clear();

    auto pushUnique = [](std::vector<std::string>& v, const std::string& s) {
        if (std::find(v.begin(), v.end(), s) == v.end()) {
            v.push_back(s);
        }
    };

    std::string line;
    while (std::getline(file, line)) {
        line = InpParserStructures::trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto tokens = InpParserStructures::split(line, ':');
        if (tokens.size() != 2) continue;

        std::string key = InpParserStructures::trim(tokens[0]);
        std::string value = InpParserStructures::trim(tokens[1]);

        pushUnique(basic_info.inp_keys_seen, key);

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

    auto& basic_info = param.get_basic_info();
    basic_info.inp_keys_seen.clear();
    basic_info.inp_keys_handled.clear();
    basic_info.inp_keys_unknown.clear();
    basic_info.inp_keys_ignored.clear();

    auto pushUnique = [](std::vector<std::string>& v, const std::string& s) {
        if (std::find(v.begin(), v.end(), s) == v.end()) {
            v.push_back(s);
        }
    };

    while (std::getline(iss, line)) {
        line = InpParserStructures::trim(line);
        if (line.empty() || line[0] == '#') continue;

        auto tokens = InpParserStructures::split(line, ':');
        if (tokens.size() != 2) continue;

        std::string key = InpParserStructures::trim(tokens[0]);
        std::string value = InpParserStructures::trim(tokens[1]);

        pushUnique(basic_info.inp_keys_seen, key);

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
    bool handled = false;

    auto pushUnique = [](std::vector<std::string>& v, const std::string& s) {
        if (std::find(v.begin(), v.end(), s) == v.end()) {
            v.push_back(s);
        }
    };

    // File paths
    if (key == "par") {
        handled = true;
        file_info.par_files.push_back(value);
    } else if (key == "fragmqtr") {
        handled = true;
        // Legacy gcmc_gpu key: additional MQTR input files
        file_info.fragment_mqtr_files.push_back(value);
        // MQTR functionality is not implemented in gcmc_cpu yet.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "fragitp") {
        handled = true;
        file_info.fragment_top_files.push_back(value);
    } else if (key == "atomtypes") {
        handled = true;
        file_info.atomtype_file = value;
        // Legacy compatibility key; currently not used by gcmc_cpu.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "monomerdir") {
        handled = true;
        file_info.monomer_dir = value;
    } else if (key == "top") {
        handled = true;
        file_info.topology_file = value;
    } else if (key == "pdb") {
        handled = true;
        file_info.input_pdb_file = value;
    } else if (key == "protitp") {
        handled = true;
        file_info.protein_top_files.push_back(value);
        // Legacy compatibility key; currently not used by gcmc_cpu.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "op_top") {
        handled = true;
        file_info.output_top_file = value;
    } else if (key == "op_pdb") {
        handled = true;
        file_info.output_pdb_file = value;
    } else if (key == "conc_norm") {
        handled = true;
        file_info.conc_norm = value;
        // Legacy compatibility key; currently not used by gcmc_cpu.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "conc_region") {
        handled = true;
        file_info.conc_region = value;
        // Legacy compatibility key; currently not used by gcmc_cpu.
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "inp_units" || key == "units") {
        handled = true;
        // Store raw unit system string for later conversion in InpParserGCMC::enhance_param
        std::string v = value;
        std::transform(v.begin(), v.end(), v.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        basic_info.inp_units = v;
        basic_info.inp_units_explicit = true;
    } else if (key == "version") {
        handled = true;
        // Keep raw version string; InpParserGCMC may use this as a hint for legacy gcmc_gpu unit mode.
        basic_info.version = value;
    } else if (key == "random_seed" || key == "seed") {
        handled = true;
        // Prefer compatibility with gcmc_gpu's "random_seed" key.
        // Use 0 as "auto" (matches existing behavior for unsigned random_seed).
        try {
            const int s = std::stoi(value);
            basic_info.random_seed = s > 0 ? static_cast<unsigned int>(s) : 0u;
        } catch (const std::exception&) {
            // Ignore invalid seed values
        }
    }
    // Space parameters
    else if (key == "grid_dx") {
        handled = true;
        // Raw value; normalized to internal nm in InpParserGCMC::enhance_param.
        space_info.grid_spacing = std::stof(value);
    } else if (key == "box_size" || key == "box") {
        handled = true;
        space_info.box_size = InpParserStructures::parse_float_array(value);
        // Raw value; normalized to internal nm in InpParserGCMC::enhance_param.
        // Calculate volume in the same raw units; it will be normalized later.
        space_info.volume = space_info.box_size[0] * space_info.box_size[1] * space_info.box_size[2];
    } else if (key == "cutoff") {
        handled = true;
        // Raw value; normalized to internal nm in InpParserGCMC::enhance_param.
        space_info.cutoff = std::stof(value);
        space_info.cutoff_explicit = true;
    } else if (key == "gc_center") {
        handled = true;
        space_info.gc_center = InpParserStructures::parse_float_array(value);
        // Raw value; normalized to internal nm in InpParserGCMC::enhance_param.
    } else if (key == "sys_center") {
        handled = true;
        space_info.sys_center = InpParserStructures::parse_float_array(value);
        // Raw value; normalized to internal nm in InpParserGCMC::enhance_param.
    }
    // Fragment parameters - support both single and multiple entries
    else if (key == "fragname") {
        handled = true;
        // Support both single fragment and comma-separated list
        auto names = InpParserStructures::parse_string_vector(value);
        for (const auto& name : names) {
            file_info.fragment_names.push_back(name);
        }
    } else if (key == "fragconc") {
        handled = true;
        // For single-component systems, override previous values
        // For multi-component, accumulate values
        auto concs = InpParserStructures::parse_float_vector(value);
        if (file_info.fragment_names.size() <= 1 && !fragment_info.conc_list.empty()) {
            // Single component mode - override
            fragment_info.conc_list.clear();
        }
        for (float conc : concs) {
            fragment_info.conc_list.push_back(conc);
        }
    } else if (key == "fragmuex") {
        handled = true;
        // For single-component systems, override previous values
        // For multi-component, accumulate values
        auto muexs = InpParserStructures::parse_float_vector(value);
        if (file_info.fragment_names.size() <= 1 && !fragment_info.muex_list.empty()) {
            // Single component mode - override
            fragment_info.muex_list.clear();
        }
        for (float muex : muexs) {
            fragment_info.muex_list.push_back(muex);
        }
    }
    // MC parameters
    else if (key == "nprint") {
        handled = true;
        mc_info.print_freq = std::stoi(value);
    } else if (key == "nsave") {
        handled = true;
        mc_info.save_freq = std::stoi(value);
    } else if (key == "mcsteps") {
        handled = true;
        mc_info.mc_steps = std::stoi(value);
    } else if (key == "moves_per_step" || key == "movesPerStep") {
        handled = true;
        mc_info.moves_per_step = std::stoi(value);
    } else if (key == "temperature") {
        handled = true;
        mc_info.temperature = std::stof(value);
        // Calculate beta from temperature (beta = 1/(kB*T) in kJ/mol units)
        mc_info.beta = 1.0f / (mc_info.BOLTZMANN * mc_info.temperature);
    } else if (key == "eqsteps") {
        handled = true;
        // Store equilibration steps if needed
        // Currently not used in MCParams, but parsed for compatibility
        pushUnique(basic_info.inp_keys_ignored, key);
    }
    // Bias parameters
    else if (key == "use_cavity_bias") {
        handled = true;
        bias_info.use_cavity_bias = (value == "yes");
    } else if (key == "use_conf_bias") {
        handled = true;
        bias_info.use_conf_bias = (value == "yes");
    }
    // Basic info parameters
    else if (key == "initcycle") {
        handled = true;
        basic_info.init_cycle = (value == "yes");
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "conserve_frags") {
        handled = true;
        basic_info.conserve_fragments = (value == "yes");
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "map_generation") {
        handled = true;
        file_info.generate_maps = (value == "yes");
        pushUnique(basic_info.inp_keys_ignored, key);
    } else if (key == "map_filename_prefix") {
        handled = true;
        file_info.map_prefix = value;
        pushUnique(basic_info.inp_keys_ignored, key);
    }

    if (handled) {
        pushUnique(basic_info.inp_keys_handled, key);
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
    // Both fragconc and fragmuex are optional
    if (!fragment_info.conc_list.empty() &&
        file_info.fragment_names.size() != fragment_info.conc_list.size()) {
        throw std::runtime_error("Inconsistent fragment parameters: fragname and fragconc sizes don't match");
    }
    if (!fragment_info.muex_list.empty() &&
        file_info.fragment_names.size() != fragment_info.muex_list.size()) {
        throw std::runtime_error("Inconsistent fragment parameters: fragname and fragmuex sizes don't match");
    }

    // At least one of fragconc or fragmuex must be specified
    if (fragment_info.conc_list.empty() && fragment_info.muex_list.empty()) {
        throw std::runtime_error("Must specify either fragconc or fragmuex (or both) for GCMC");
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
