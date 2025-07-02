// src/io/topology/psfParserMain.cpp

#include "psfParserMain.hpp"
#include "psfParserStringUtils.hpp"
#include "psfParserSections.hpp"
#include "psfParserSectionsConnectivity.hpp"
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <set>
#include <filesystem>

namespace pygcmc {
namespace io {

model::Topology PSFParser::parse_file(const std::string& filename) {
    model::Topology topology;
    PSFParser parser;
    if (!parser.parse_to_topology(filename, topology)) {
        throw std::runtime_error("Failed to parse PSF file: " + filename);
    }
    return topology;
}

model::Topology PSFParser::parse_string(const std::string& psf_str) {
    return PSFParserStringUtils::parse_string(psf_str);
}

// Helper function to read file lines
bool readFileToLines(const std::string& filename, std::vector<std::string>& lines) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        lines.push_back(line);
    }

    return true;
}

std::string PSFParser::trim(const std::string& str) {
    return PSFParserStringUtils::trim(str);
}

bool PSFParser::parse_to_topology(const std::string& filename, model::Topology& topology) {
    // Read all lines at once
    std::vector<std::string> lines;
    if (!readFileToLines(filename, lines)) {
        return false;
    }
    
    // Check for extended Drude format in header
    bool is_extended_format = false;
    bool is_drude_format = false;
    if (!lines.empty()) {
        std::string header = trim(lines[0]);
        if (header.find("PSF") != std::string::npos) {
            is_extended_format = header.find("EXT") != std::string::npos;
            is_drude_format = header.find("DRUDE") != std::string::npos;
            
            if (is_drude_format && is_extended_format) {
                std::cout << "Detected extended Drude PSF format" << std::endl;
            }
        }
    }

    // First, collect all sections
    std::unordered_map<std::string, std::vector<std::string>> sections;
    std::string current_section;
    std::vector<std::string> current_lines;

    for (const auto& line : lines) {
        std::string trimmed = trim(line);
        if (trimmed.empty()) continue;

        // Check if this is a section header
        if (trimmed.find('!') != std::string::npos) {
            // Save previous section if any
            if (!current_section.empty() && !current_lines.empty()) {
                sections[current_section] = current_lines;
            }

            // Determine new section type
            if (trimmed.find("!NATOM") != std::string::npos) {
                current_section = "NATOM";
            } else if (trimmed.find("!NBOND") != std::string::npos) {
                current_section = "NBOND";
            } else if (trimmed.find("!NTHETA") != std::string::npos) {
                current_section = "NTHETA";
            } else if (trimmed.find("!NPHI") != std::string::npos) {
                current_section = "NPHI";
            } else if (trimmed.find("!NIMPHI") != std::string::npos) {
                current_section = "NIMPHI";
            } else if (trimmed.find("!NCRTERM") != std::string::npos ||
                      trimmed.find("!NCMAP") != std::string::npos) {
                current_section = "CMAP";
            } else if (trimmed.find("!NGRP") != std::string::npos) {
                current_section = "NGRP";
            } else if (trimmed.find("!NDON") != std::string::npos) {
                current_section = "NDON";
            } else if (trimmed.find("!NACC") != std::string::npos) {
                current_section = "NACC";
            } else {
                current_section = "";  // Unknown section
            }
            current_lines.clear();
            current_lines.push_back(trimmed);  // Include the header line
        } else if (!current_section.empty()) {
            // Add line to current section
            current_lines.push_back(trimmed);
        }
    }

    // Add the last section if any
    if (!current_section.empty() && !current_lines.empty()) {
        sections[current_section] = current_lines;
    }

    // Now process sections in the required order
    // NATOM must be first
    if (!sections.count("NATOM")) {
        std::cerr << "Missing required NATOM section" << std::endl;
        return false;
    }

    if (!PSFParserSections::parse_atoms_from_lines(sections["NATOM"], topology, is_extended_format, is_drude_format)) {
        std::cerr << "Failed to parse NATOM section" << std::endl;
        return false;
    }

    // Process optional sections in order
    const std::vector<std::string> optional_sections = {
        "NBOND", "NTHETA", "NPHI", "NIMPHI", "CMAP", "NGRP", "NDON", "NACC"
    };

    for (const auto& section : optional_sections) {
        if (sections.count(section)) {
            bool success = false;

            if (section == "NBOND") {
                success = PSFParserSectionsConnectivity::parse_bonds_from_lines(sections[section], topology);
            } else if (section == "NTHETA") {
                success = PSFParserSectionsConnectivity::parse_angles_from_lines(sections[section], topology);
            } else if (section == "NPHI") {
                success = PSFParserSectionsConnectivity::parse_dihedrals_from_lines(sections[section], topology);
            } else if (section == "NIMPHI") {
                success = PSFParserSections::parse_impropers_from_lines(sections[section], topology);
            } else if (section == "CMAP") {
                success = PSFParserSections::parse_cmap_from_lines(sections[section], topology);
            } else if (section == "NGRP") {
                success = PSFParserSections::parse_groups_from_lines(sections[section], topology);
            } else if (section == "NDON") {
                success = PSFParserSections::parse_donors_from_lines(sections[section], topology);
            } else if (section == "NACC") {
                success = PSFParserSections::parse_acceptors_from_lines(sections[section], topology);
            }

            if (!success) {
                std::cerr << "Failed to parse " << section << " section" << std::endl;
                return false;
            }
        }
    }

    return true;
}

} // namespace io
} // namespace pygcmc