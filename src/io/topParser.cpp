// src/io/topParser.cpp

#include "topParser.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <filesystem>
#include <stdexcept>

namespace pygcmc {
namespace io {

model::Topology TopParser::parse(const std::string& filename) {
    model::Topology topology;
    processed_files_.clear();
    current_molecule_type_.clear();
    molecule_order_.clear();
    
    // First process includes
    process_includes(filename, topology);
    
    // Parse all sections in order
    auto defaults_lines = read_section(filename, "defaults");
    auto atomtypes_lines = read_section(filename, "atomtypes");
    auto moleculetype_lines = read_section(filename, "moleculetype");
    auto atoms_lines = read_section(filename, "atoms");
    auto bonds_lines = read_section(filename, "bonds");
    auto angles_lines = read_section(filename, "angles");
    auto dihedrals_lines = read_section(filename, "dihedrals");
    auto impropers_lines = read_section(filename, "impropers");
    auto pairs_lines = read_section(filename, "pairs");
    auto exclusions_lines = read_section(filename, "exclusions");
    auto cmap_lines = read_section(filename, "cmap");
    auto system_lines = read_section(filename, "system");
    auto molecules_lines = read_section(filename, "molecules");
    
    // Parse each section if it exists
    if (!defaults_lines.empty()) parse_defaults_section(defaults_lines, topology);
    if (!atomtypes_lines.empty()) parse_atomtypes_section(atomtypes_lines, topology);
    if (!moleculetype_lines.empty()) parse_moleculetype_section(moleculetype_lines, topology);
    if (!atoms_lines.empty()) parse_atoms_section(atoms_lines, topology);
    if (!bonds_lines.empty()) parse_bonds_section(bonds_lines, topology);
    if (!angles_lines.empty()) parse_angles_section(angles_lines, topology);
    if (!dihedrals_lines.empty()) parse_dihedrals_section(dihedrals_lines, topology);
    if (!impropers_lines.empty()) parse_impropers_section(impropers_lines, topology);
    if (!pairs_lines.empty()) parse_pairs_section(pairs_lines, topology);
    if (!exclusions_lines.empty()) parse_exclusions_section(exclusions_lines, topology);
    if (!cmap_lines.empty()) parse_cmap_section(cmap_lines, topology);
    if (!system_lines.empty()) parse_system_section(system_lines, topology);
    if (!molecules_lines.empty()) parse_molecules_section(molecules_lines, topology);
    
    // Validate the topology after parsing
    if (!validate_topology(topology)) {
        report_error("Topology validation failed", true);
    }
    
    return topology;
}

bool TopParser::process_includes(const std::string& filename, model::Topology& topology) {
    // Check if we've already processed this file
    if (processed_files_.find(filename) != processed_files_.end()) {
        return true;  // Skip already processed files
    }
    
    processed_files_.insert(filename);
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        report_error("Could not open file for include processing: " + filename, false);
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        trim(line);
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == ';') {
            continue;
        }

        // Handle preprocessor directives
        if (line[0] == '#') {
            if (!handle_preprocessor_line(line)) {
                continue;
            }
        }
        
        // Skip lines if we're in a false #ifdef block
        if (!should_process_line()) {
            continue;
        }
    }

    return true;
}

bool TopParser::handle_preprocessor_line(const std::string& line) {
    std::istringstream iss(line);
    std::string directive;
    iss >> directive;
    
    if (directive == "#include") {
        std::string include_path;
        iss >> include_path;
        // Remove quotes if present
        if (include_path.front() == '"' && include_path.back() == '"') {
            include_path = include_path.substr(1, include_path.length() - 2);
        }
        
        std::string resolved_path = resolve_include_path(include_path, "");
        if (!resolved_path.empty()) {
            return process_includes(resolved_path, topology);
        }
        return false;
    }
    else if (directive == "#ifdef") {
        std::string condition;
        iss >> condition;
        preproc_state_.in_ifdef = true;
        preproc_state_.ifdef_depth++;
        preproc_state_.ifdef_condition_met = evaluate_ifdef_condition(condition);
        preproc_state_.current_ifdef = condition;
        return true;
    }
    else if (directive == "#ifndef") {
        std::string condition;
        iss >> condition;
        preproc_state_.in_ifdef = true;
        preproc_state_.ifdef_depth++;
        preproc_state_.ifdef_condition_met = !evaluate_ifdef_condition(condition);
        preproc_state_.current_ifdef = condition;
        return true;
    }
    else if (directive == "#else") {
        if (preproc_state_.in_ifdef) {
            preproc_state_.in_else = true;
            preproc_state_.ifdef_condition_met = !preproc_state_.ifdef_condition_met;
        }
        return true;
    }
    else if (directive == "#endif") {
        if (preproc_state_.ifdef_depth > 0) {
            preproc_state_.ifdef_depth--;
            if (preproc_state_.ifdef_depth == 0) {
                preproc_state_.in_ifdef = false;
                preproc_state_.in_else = false;
                preproc_state_.ifdef_condition_met = true;
                preproc_state_.current_ifdef.clear();
            }
        }
        return true;
    }
    
    return true;
}

bool TopParser::evaluate_ifdef_condition(const std::string& condition) {
    return preprocessor_defines_.find(condition) != preprocessor_defines_.end();
}

bool TopParser::should_process_line() const {
    if (!preproc_state_.in_ifdef) {
        return true;
    }
    return preproc_state_.ifdef_condition_met;
}

void TopParser::report_error(const std::string& message, bool critical) {
    if (critical && strict_mode_) {
        throw std::runtime_error(message);
    }
    std::cerr << (critical ? "Error: " : "Warning: ") << message << std::endl;
}

bool TopParser::validate_topology(const model::Topology& topology) {
    // Check molecule consistency
    if (!check_molecule_consistency(topology)) {
        return false;
    }
    
    // Add more validation as needed
    return true;
}

bool TopParser::check_molecule_consistency(const model::Topology& topology) {
    // Check if the number of atoms matches the molecule definitions
    size_t total_atoms = 0;
    for (const auto& mol_count : molecule_order_) {
        const auto& mol_type = topology.get_moleculetypes().find(mol_count.first);
        if (mol_type == topology.get_moleculetypes().end()) {
            report_error("Molecule type " + mol_count.first + " not found in topology", true);
            return false;
        }
        total_atoms += mol_type->second.atoms.size() * mol_count.second;
    }
    
    if (total_atoms != topology.get_num_atoms()) {
        report_error("Total atom count mismatch: expected " + std::to_string(total_atoms) + 
                    " but got " + std::to_string(topology.get_num_atoms()), true);
        return false;
    }
    
    return true;
}

std::string TopParser::resolve_include_path(const std::string& include_path, const std::string& parent_file) {
    namespace fs = std::filesystem;
    
    // First try relative to the parent file
    fs::path parent_dir = fs::path(parent_file).parent_path();
    fs::path resolved = parent_dir / include_path;
    
    if (fs::exists(resolved)) {
        return resolved.string();
    }
    
    // If not found, try some common locations
    std::vector<fs::path> search_paths = {
        fs::current_path(),
        fs::current_path() / "include",
        fs::current_path() / "forcefield",
        parent_dir / "include",
        parent_dir / "forcefield"
    };
    
    for (const auto& path : search_paths) {
        resolved = path / include_path;
        if (fs::exists(resolved)) {
            return resolved.string();
        }
    }
    
    std::cerr << "Warning: Include file not found: " << include_path << std::endl;
    return "";
}

std::vector<std::string> TopParser::read_section(const std::string& filename, const std::string& section_name) {
    std::vector<std::string> lines;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return lines;
    }

    std::string line;
    bool in_section = false;
    std::string section_header = "[ " + section_name + " ]";

    while (std::getline(file, line)) {
        trim(line);
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == ';' || line[0] == '#') {
            continue;
        }

        // Check for section start
        if (line == section_header) {
            in_section = true;
            continue;
        }
        
        // Check for section end (new section starts)
        if (in_section && line[0] == '[') {
            break;
        }

        // Store lines within the section
        if (in_section) {
            lines.push_back(line);
        }
    }

    return lines;
}

bool TopParser::parse_defaults_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 5) continue;

        try {
            model::TopologyDefaults defaults;
            defaults.nbfunc = std::stoi(tokens[0]);
            defaults.combination_rule = std::stoi(tokens[1]);
            defaults.gen_pairs = (tokens[2] == "yes" || tokens[2] == "YES");
            defaults.fudgeLJ = std::stod(tokens[3]);
            defaults.fudgeQQ = std::stod(tokens[4]);
            
            topology.set_defaults(defaults);
            return true;  // Only need to process the first valid line
        } catch (const std::exception& e) {
            std::cerr << "Error parsing defaults line: " << line << std::endl;
            continue;
        }
    }
    return false;
}

bool TopParser::parse_atomtypes_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 7) continue;

        try {
            model::AtomType atomtype;
            atomtype.name = tokens[0];
            atomtype.atomic_number = tokens.size() > 1 ? std::stoi(tokens[1]) : 0;
            atomtype.mass = std::stod(tokens[2]);
            atomtype.charge = std::stod(tokens[3]);
            atomtype.ptype = tokens[4][0];  // Take first character
            atomtype.sigma = std::stod(tokens[5]);
            atomtype.epsilon = std::stod(tokens[6]);
            
            topology.add_atomtype(atomtype);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing atomtype line: " << line << std::endl;
            continue;
        }
    }
    return true;
}

bool TopParser::parse_moleculetype_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 2) continue;

        try {
            current_molecule_type_ = tokens[0];
            current_molecule_nrexcl_ = std::stoi(tokens[1]);
            
            model::MoleculeType moltype;
            moltype.name = current_molecule_type_;
            moltype.nrexcl = current_molecule_nrexcl_;
            
            topology.add_moleculetype(moltype);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing moleculetype line: " << line << std::endl;
            continue;
        }
    }
    return true;
}

bool TopParser::parse_atoms_section(const std::vector<std::string>& lines, model::Topology& topology) {
    std::map<int, std::string> segment_map;  // residue_number -> segment_name
    
    // Determine segment type based on current molecule type
    std::string current_segment = molecule_to_segment_type_[current_molecule_type_];
    if (current_segment.empty()) {
        current_segment = "PROT";  // Default segment name
    }
    
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 8) continue;  // Need at least 8 columns

        try {
            // Reuse existing TopologyAtom structure through add_atom
            int atom_id = std::stoi(tokens[0]) - 1;  // Convert to 0-based indexing
            std::string atom_type = tokens[1];
            int residue_number = std::stoi(tokens[2]);
            std::string residue_name = tokens[3];
            std::string atom_name = tokens[4];
            double charge = std::stod(tokens[6]);
            double mass = std::stod(tokens[7]);

            // Add atom to topology using existing PSF data structures
            topology.add_atom(atom_name, atom_type, charge, mass,
                            residue_name, residue_number, current_segment);

            // Store B-state parameters if present (we'll need to extend TopologyAtom for this)
            if (tokens.size() >= 10) {
                // TODO: Add B-state support to TopologyAtom
                // For now, we'll store this information in the molecule type
                if (auto mol_it = topology.get_moleculetypes().find(current_molecule_type_);
                    mol_it != topology.get_moleculetypes().end()) {
                    // Store B-state info in molecule type
                    mol_it->second.atoms.push_back(atom_id);
                }
            }

        } catch (const std::exception& e) {
            report_error("Error parsing atom line: " + line + "\n" + e.what(), false);
            continue;
        }
    }
    return true;
}

bool TopParser::parse_bonds_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 3) continue;  // Need at least ai, aj, funct

        try {
            int atom1 = std::stoi(tokens[0]) - 1;
            int atom2 = std::stoi(tokens[1]) - 1;
            int funct = std::stoi(tokens[2]);
            
            double length = 0.0;
            double force_constant = 0.0;
            
            if (tokens.size() >= 5) {
                length = std::stod(tokens[3]);
                force_constant = std::stod(tokens[4]);
            }
            
            // Pass function_type directly to add_bond
            topology.add_bond(atom1, atom2, length, force_constant, funct);

        } catch (const std::exception& e) {
            report_error("Error parsing bond line: " + line, false);
            continue;
        }
    }
    return true;
}

bool TopParser::parse_angles_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 4) continue;  // Need ai, aj, ak, funct

        try {
            int atom1 = std::stoi(tokens[0]) - 1;
            int atom2 = std::stoi(tokens[1]) - 1;
            int atom3 = std::stoi(tokens[2]) - 1;
            int funct = std::stoi(tokens[3]);
            
            double angle = 0.0;
            double force_constant = 0.0;
            
            if (tokens.size() >= 6) {
                angle = std::stod(tokens[4]);
                force_constant = std::stod(tokens[5]);
            }
            
            // Pass function_type directly to add_angle
            topology.add_angle(atom1, atom2, atom3, angle, force_constant, funct);

        } catch (const std::exception& e) {
            report_error("Error parsing angle line: " + line, false);
            continue;
        }
    }
    return true;
}

bool TopParser::parse_dihedrals_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 5) continue;  // Need ai, aj, ak, al, funct

        try {
            int atom1 = std::stoi(tokens[0]) - 1;
            int atom2 = std::stoi(tokens[1]) - 1;
            int atom3 = std::stoi(tokens[2]) - 1;
            int atom4 = std::stoi(tokens[3]) - 1;
            int funct = std::stoi(tokens[4]);
            
            int multiplicity = 1;
            double angle = 0.0;
            double force_constant = 0.0;
            
            if (tokens.size() >= 8) {
                multiplicity = std::stoi(tokens[5]);
                angle = std::stod(tokens[6]);
                force_constant = std::stod(tokens[7]);
            }
            
            // Pass function_type directly to add_dihedral
            topology.add_dihedral(atom1, atom2, atom3, atom4, multiplicity,
                                angle, force_constant, false, funct);  // not improper

        } catch (const std::exception& e) {
            report_error("Error parsing dihedral line: " + line, false);
            continue;
        }
    }
    return true;
}

bool TopParser::parse_impropers_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 5) continue;  // Need ai, aj, ak, al, funct

        try {
            int atom1 = std::stoi(tokens[0]) - 1;
            int atom2 = std::stoi(tokens[1]) - 1;
            int atom3 = std::stoi(tokens[2]) - 1;
            int atom4 = std::stoi(tokens[3]) - 1;
            int funct = std::stoi(tokens[4]);
            
            double angle = 0.0;
            double force_constant = 0.0;
            
            if (tokens.size() >= 7) {
                angle = std::stod(tokens[5]);
                force_constant = std::stod(tokens[6]);
            }
            
            // Pass function_type directly to add_dihedral with improper=true
            topology.add_dihedral(atom1, atom2, atom3, atom4, 0,
                                angle, force_constant, true, funct);  // improper = true

        } catch (const std::exception& e) {
            report_error("Error parsing improper line: " + line, false);
            continue;
        }
    }
    return true;
}

bool TopParser::parse_pairs_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 2) continue;

        try {
            int atom1 = std::stoi(tokens[0]) - 1;  // Convert to 0-based indexing
            int atom2 = std::stoi(tokens[1]) - 1;
            
            double c6 = 0.0;
            double c12 = 0.0;
            
            if (tokens.size() >= 4) {
                c6 = std::stod(tokens[3]);
                if (tokens.size() >= 5) {
                    c12 = std::stod(tokens[4]);
                }
            }
            
            topology.add_nonbonded_exclusion(atom1, atom2);  // Pairs are also treated as exclusions

        } catch (const std::exception& e) {
            std::cerr << "Error parsing pair line: " << line << std::endl;
            continue;
        }
    }
    return true;
}

bool TopParser::parse_exclusions_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 2) continue;

        try {
            int atom1 = std::stoi(tokens[0]) - 1;  // Convert to 0-based indexing
            
            // An exclusion line can contain multiple atoms to be excluded from atom1
            for (size_t i = 1; i < tokens.size(); ++i) {
                int atom2 = std::stoi(tokens[i]) - 1;
                topology.add_nonbonded_exclusion(atom1, atom2);
            }

        } catch (const std::exception& e) {
            std::cerr << "Error parsing exclusion line: " << line << std::endl;
            continue;
        }
    }
    return true;
}

bool TopParser::parse_cmap_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 9) continue;  // Need 8 atoms + funct

        try {
            model::CMAPTerm cmap;
            for (int i = 0; i < 8; ++i) {
                cmap.atoms[i] = std::stoi(tokens[i]) - 1;
            }
            // Note: actual CMAP grid data would typically be in a separate file
            // Here we just record the atom indices
            topology.add_cmap(cmap);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing CMAP line: " << line << std::endl;
            continue;
        }
    }
    return true;
}

bool TopParser::parse_system_section(const std::vector<std::string>& lines, model::Topology& topology) {
    // System section is typically just a title, no need to parse
    return true;
}

bool TopParser::parse_molecules_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 2) continue;

        try {
            std::string mol_name = tokens[0];
            int count = std::stoi(tokens[1]);
            
            // Store molecule order for validation
            molecule_order_.emplace_back(mol_name, count);
            
            // Determine segment type based on molecule name
            std::string segment_type = "PROT";  // Default
            if (mol_name.find("SOL") != std::string::npos || 
                mol_name.find("WAT") != std::string::npos) {
                segment_type = "SOL";
            }
            else if (mol_name.find("ION") != std::string::npos) {
                segment_type = "ION";
            }
            else if (mol_name.find("LIG") != std::string::npos) {
                segment_type = "LIG";
            }
            molecule_to_segment_type_[mol_name] = segment_type;
            
            topology.add_molecule(mol_name, count);
        } catch (const std::exception& e) {
            report_error("Error parsing molecules line: " + line, false);
            continue;
        }
    }
    return true;
}

void TopParser::trim(std::string& str) {
    // Trim leading spaces
    str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
    
    // Trim trailing spaces
    str.erase(std::find_if(str.rbegin(), str.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), str.end());
}

std::vector<std::string> TopParser::split(const std::string& str) {
    std::vector<std::string> tokens;
    std::istringstream iss(str);
    std::string token;
    
    while (iss >> token) {
        if (token[0] == ';') break;  // Stop at comments
        tokens.push_back(token);
    }
    
    return tokens;
}

} // namespace io
} // namespace pygcmc


