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

bool TopParser::parse_to_topology(const std::string& filename, model::Topology& topology) {
    processed_files_.clear();
    current_molecule_type_.clear();
    molecule_order_.clear();
    
    // First process includes
    if (!process_includes(filename, topology)) {
        return false;
    }
    
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
    
    // Validate that at least one required section exists
    if (moleculetype_lines.empty() && atoms_lines.empty()) {
        std::cerr << "Invalid topology file: no moleculetype or atoms sections found" << std::endl;
        return false;
    }
    
    // Parse each section if it exists
    // Skip defaults and atomtypes for now as they are force field parameters
    if (!moleculetype_lines.empty() && !parse_moleculetype_section(moleculetype_lines, topology)) return false;
    if (!atoms_lines.empty() && !parse_atoms_section(atoms_lines, topology)) return false;
    if (!bonds_lines.empty() && !parse_bonds_section(bonds_lines, topology)) return false;
    if (!angles_lines.empty() && !parse_angles_section(angles_lines, topology)) return false;
    if (!dihedrals_lines.empty() && !parse_dihedrals_section(dihedrals_lines, topology)) return false;
    if (!impropers_lines.empty() && !parse_impropers_section(impropers_lines, topology)) return false;
    if (!pairs_lines.empty() && !parse_pairs_section(pairs_lines, topology)) return false;
    if (!exclusions_lines.empty() && !parse_exclusions_section(exclusions_lines, topology)) return false;
    if (!cmap_lines.empty() && !parse_cmap_section(cmap_lines, topology)) return false;
    if (!system_lines.empty() && !parse_system_section(system_lines, topology)) return false;
    if (!molecules_lines.empty() && !parse_molecules_section(molecules_lines, topology)) return false;
    
    return true;
}

bool TopParser::process_includes(const std::string& filename, model::Topology& topology) {
    // Check if we've already processed this file
    if (processed_files_.find(filename) != processed_files_.end()) {
        return true;  // Skip already processed files
    }
    
    processed_files_.insert(filename);
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Could not open file for include processing: " << filename << std::endl;
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
                
                std::string resolved_path = resolve_include_path(include_path, filename);
                if (!resolved_path.empty()) {
                    if (!process_includes(resolved_path, topology)) {
                        return false;
                    }
                }
            }
            // Skip other preprocessor directives for now
            continue;
        }
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

std::string TopParser::trim(std::string& str) {
    // Trim leading spaces
    str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
    
    // Trim trailing spaces
    str.erase(std::find_if(str.rbegin(), str.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), str.end());
    
    return str;
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

bool TopParser::parse_moleculetype_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 2) continue;

        try {
            current_molecule_type_ = tokens[0];
            current_molecule_nrexcl_ = std::stoi(tokens[1]);
            
            // Add as a segment since we don't have molecule types in Topology
            topology.add_segment(current_molecule_type_);
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error parsing moleculetype line: " << line << std::endl;
            return false;
        }
    }
    return true;
}

bool TopParser::parse_atoms_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 8) continue;  // Need at least 8 columns

        try {
            // Parse atom data
            std::string atom_type = tokens[1];
            int residue_number = std::stoi(tokens[2]);
            std::string residue_name = tokens[3];
            std::string atom_name = tokens[4];
            double charge = std::stod(tokens[6]);
            double mass = std::stod(tokens[7]);

            // Add atom to topology
            topology.add_atom(
                atom_name,
                atom_type,
                charge,
                mass,
                residue_name,
                residue_number,
                current_molecule_type_  // Use molecule type as segment
            );
        } catch (const std::exception& e) {
            std::cerr << "Error parsing atom line: " << line << std::endl;
            return false;
        }
    }
    return true;
}

bool TopParser::parse_bonds_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 2) continue;  // Need at least ai, aj

        try {
            int atom1 = std::stoi(tokens[0]) - 1;  // Convert to 0-based indexing
            int atom2 = std::stoi(tokens[1]) - 1;
            
            // Add bond to topology
            topology.add_bond(atom1, atom2);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing bond line: " << line << std::endl;
            return false;
        }
    }
    return true;
}

bool TopParser::parse_angles_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 3) continue;  // Need at least ai, aj, ak

        try {
            int atom1 = std::stoi(tokens[0]) - 1;  // Convert to 0-based indexing
            int atom2 = std::stoi(tokens[1]) - 1;
            int atom3 = std::stoi(tokens[2]) - 1;
            
            // Add angle to topology
            topology.add_angle(atom1, atom2, atom3);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing angle line: " << line << std::endl;
            return false;
        }
    }
    return true;
}

bool TopParser::parse_dihedrals_section(const std::vector<std::string>& lines, model::Topology& topology) {
    std::cerr << "Parsing " << lines.size() << " dihedral lines" << std::endl;
    int proper_count = 0;
    int improper_count = 0;

    for (const auto& line : lines) {
        auto tokens = split(line);
        // GROMACS dihedral lines can have 4-11 columns, but we need at least 4 atoms
        if (tokens.size() < 4) {
            std::cerr << "Skipping line with insufficient tokens: " << line << std::endl;
            continue;
        }

        try {
            int atom1 = std::stoi(tokens[0]) - 1;  // Convert to 0-based indexing
            int atom2 = std::stoi(tokens[1]) - 1;
            int atom3 = std::stoi(tokens[2]) - 1;
            int atom4 = std::stoi(tokens[3]) - 1;
            
            // Default to function type 1 (proper dihedral) if not specified
            int funcType = (tokens.size() >= 5) ? std::stoi(tokens[4]) : 1;
            
            // GROMACS function types 2 and 4 are impropers
            // Type 9 is multiple proper dihedrals
            if (funcType == 2 || funcType == 4) {
                topology.add_improper(atom1, atom2, atom3, atom4);
                improper_count++;
            } else {
                topology.add_dihedral(atom1, atom2, atom3, atom4);
                proper_count++;
                // For type 9 (multiple), we need to add additional proper dihedrals
                if (funcType == 9 && tokens.size() >= 7) {
                    int multiplicity = std::stoi(tokens[6]);
                    for (int i = 1; i < multiplicity; i++) {
                        topology.add_dihedral(atom1, atom2, atom3, atom4);
                        proper_count++;
                    }
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Error parsing dihedral line: " << line << " - " << e.what() << std::endl;
            return false;
        }
    }
    
    std::cerr << "Added " << proper_count << " proper dihedrals and " 
              << improper_count << " improper dihedrals" << std::endl;
    return true;
}

bool TopParser::parse_impropers_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 4) continue;  // Need at least ai, aj, ak, al

        try {
            int atom1 = std::stoi(tokens[0]) - 1;  // Convert to 0-based indexing
            int atom2 = std::stoi(tokens[1]) - 1;
            int atom3 = std::stoi(tokens[2]) - 1;
            int atom4 = std::stoi(tokens[3]) - 1;
            
            // Add improper to topology
            topology.add_improper(atom1, atom2, atom3, atom4);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing improper line: " << line << std::endl;
            return false;
        }
    }
    return true;
}

bool TopParser::parse_pairs_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 2) continue;  // Need at least ai, aj

        try {
            int atom1 = std::stoi(tokens[0]) - 1;  // Convert to 0-based indexing
            int atom2 = std::stoi(tokens[1]) - 1;
            
            // Add pair as exclusion
            topology.add_nonbonded_exclusion(atom1, atom2);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing pair line: " << line << std::endl;
            return false;
        }
    }
    return true;
}

bool TopParser::parse_exclusions_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 2) continue;  // Need at least ai, aj

        try {
            int atom1 = std::stoi(tokens[0]) - 1;  // Convert to 0-based indexing
            
            // Add all exclusions for this atom
            for (size_t i = 1; i < tokens.size(); ++i) {
                int atom2 = std::stoi(tokens[i]) - 1;
                topology.add_nonbonded_exclusion(atom1, atom2);
            }
        } catch (const std::exception& e) {
            std::cerr << "Error parsing exclusion line: " << line << std::endl;
            return false;
        }
    }
    return true;
}

bool TopParser::parse_cmap_section(const std::vector<std::string>& lines, model::Topology& topology) {
    for (const auto& line : lines) {
        auto tokens = split(line);
        if (tokens.size() < 8) continue;  // Need 8 atoms for CMAP

        try {
            std::array<int, 8> atoms;
            for (int i = 0; i < 8; ++i) {
                atoms[i] = std::stoi(tokens[i]) - 1;  // Convert to 0-based indexing
            }
            
            // Add CMAP to topology
            topology.add_cmap(atoms);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing CMAP line: " << line << std::endl;
            return false;
        }
    }
    return true;
}

bool TopParser::parse_system_section([[maybe_unused]] const std::vector<std::string>& lines, 
                                   [[maybe_unused]] model::Topology& topology) {
    // System section is typically just a title, we can ignore it
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
            
            // Add segment for each molecule
            topology.add_segment(mol_name);
        } catch (const std::exception& e) {
            std::cerr << "Error parsing molecules line: " << line << std::endl;
            return false;
        }
    }
    return true;
}

} // namespace io
} // namespace pygcmc


