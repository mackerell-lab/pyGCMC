// src/io/topology/PsfParserStringUtils.cpp

#include "PsfParserStringUtils.hpp"
#include "PsfParserMain.hpp"
#include <sstream>
#include <stdexcept>
#include <vector>

namespace pygcmc {
namespace io {

model::Topology PSFParserStringUtils::parse_string(const std::string& psf_str) {
    std::vector<std::string> lines;
    std::istringstream stream(psf_str);
    std::string line;
    while (std::getline(stream, line)) {
        lines.push_back(line);
    }

    model::Topology topology;
    if (!PSFParser::parse_lines_to_topology(lines, topology)) {
        throw std::runtime_error("Failed to parse PSF string");
    }
    return topology;
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
