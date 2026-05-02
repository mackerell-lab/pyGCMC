// src/io/topology/psfParserStringUtils.cpp

#include "psfParserStringUtils.hpp"
#include "psfParserMain.hpp"
#include <fstream>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace pygcmc {
namespace io {

model::Topology PSFParserStringUtils::parse_string(const std::string& psf_str) {
    // Create a temporary file to write the string to
    std::filesystem::path temp_dir = std::filesystem::temp_directory_path();
    std::filesystem::path temp_file = temp_dir / "temp_topology.psf";

    // Write the string to the temporary file
    std::ofstream out(temp_file);
    if (!out) {
        throw std::runtime_error("Failed to create temporary file for PSF parsing");
    }
    out << psf_str;
    out.close();

    try {
        // Parse the temporary file
        model::Topology topology = PSFParser::parse_file(temp_file.string());

        // Clean up
        std::filesystem::remove(temp_file);

        return topology;
    } catch (const std::exception& e) {
        // Clean up on error
        std::filesystem::remove(temp_file);
        throw;
    }
}

std::string PSFParserStringUtils::trim(const std::string& str) {
    const auto start = str.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
        return "";
    }
    const auto end = str.find_last_not_of(" \t\r\n");
    return str.substr(start, end - start + 1);
}

} // namespace io
} // namespace pygcmc
