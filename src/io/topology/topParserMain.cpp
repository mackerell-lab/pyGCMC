// src/io/topology/topParserMain.cpp

#include "topParserMain.hpp"
#include "topParserStringUtils.hpp"
#include "topParserSectionProcessor.hpp"
#include "topParserMoleculeBuilder.hpp"
#include <stdexcept>

namespace pygcmc {
namespace io {

model::Topology TOPParser::parse_file(const std::string& filename) {
    model::Topology topology;
    TOPParser parser;
    if (!parser.parse_to_topology(filename, topology)) {
        throw std::runtime_error("Failed to parse topology file: " + filename);
    }
    return topology;
}

model::Topology TOPParser::parse_string(const std::string& top_str) {
    return TopParserStringUtils::parse_string(top_str);
}

bool TOPParser::parse_to_topology(const std::string& filename, model::Topology& topology) {
    TopParserUtilities::debug_print("\n=== Starting topology parsing of ", filename, " ===\n");

    // Clear any previous state
    molecule_order_.clear();
    processed_files_.clear();
    current_molecule_type_.clear();
    current_molecule_nrexcl_ = 0;

    // Collect all lines with preprocessor handling
    std::vector<LineInfo> all_lines;
    PreprocessorState pp_state;
    if (!TopParserPreprocessor::collect_all_lines(filename, all_lines, pp_state, true, processed_files_)) {
        return false;
    }

    if (all_lines.empty()) {
        TopParserUtilities::debug_print("Error: No valid content found in topology file\n");
        return false;
    }

    // Process sections
    SectionProcessorResult result = TopParserSectionProcessor::process_sections(all_lines);

    // Return false if no valid sections were found
    if (!result.found_valid_section) {
        TopParserUtilities::debug_print("Error: No valid topology sections found in file\n");
        return false;
    }

    // Store molecule order for later use
    molecule_order_ = result.molecule_order;

    // Build molecules and add to topology
    return TopParserMoleculeBuilder::build_molecules(topology, result);
}

} // namespace io
} // namespace pygcmc
