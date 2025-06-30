// src/io/topology/topParserStringUtils.cpp

#include "topParserStringUtils.hpp"
#include "topParserMain.hpp"
#include "topParserUtilities.hpp"
#include <fstream>
#include <filesystem>
#include <stdexcept>

namespace pygcmc {
namespace io {

model::Topology TopParserStringUtils::parse_string(const std::string& top_str) {
    // Handle empty string case - first trim whitespace
    std::string trimmed_str = top_str;
    TopParserUtilities::trim(trimmed_str);
    if (trimmed_str.empty()) {
        return model::Topology();
    }
    
    // Create a temporary file to write the string to
    std::filesystem::path temp_dir = std::filesystem::temp_directory_path();
    std::filesystem::path temp_file = temp_dir / "temp_topology.top";
    
    // Write the string to the temporary file
    std::ofstream out(temp_file);
    if (!out) {
        throw std::runtime_error("Failed to create temporary file for topology parsing");
    }
    out << top_str;  // Write original string to preserve formatting
    out.close();
    
    try {
        // Parse the temporary file
        model::Topology topology = TOPParser::parse_file(temp_file.string());
        
        // Clean up
        std::filesystem::remove(temp_file);
        
        return topology;
    } catch (const std::exception& e) {
        // Clean up on error
        std::filesystem::remove(temp_file);
        throw;
    }
}

} // namespace io
} // namespace pygcmc