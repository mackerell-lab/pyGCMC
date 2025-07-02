// src/io/structure/PdbParserMain.cpp

#include "PdbParserMain.hpp"
#include "PdbParserStructureRecords.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace pygcmc {
namespace io {
namespace structure {

model::Structure PdbParserMain::parse_file(const std::string& filename) {
    model::Structure structure;
    if (!parse_to_structure(filename, structure)) {
        throw std::runtime_error("Failed to parse PDB file: " + filename);
    }
    return structure;
}

model::Structure PdbParserMain::parse_string(const std::string& pdbStr) {
    model::Structure structure;
    if (!parse_string_to_structure(pdbStr, structure)) {
        throw std::runtime_error("Failed to parse PDB string");
    }
    return structure;
}

bool PdbParserMain::parse_to_structure(const std::string& filename, model::Structure& structure) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }
    
    return parse_common(file, structure);
}

bool PdbParserMain::parse_string_to_structure(const std::string& pdbStr, model::Structure& structure) {
    std::istringstream stream(pdbStr);
    return parse_common(stream, structure);
}

bool PdbParserMain::parse_common(std::istream& input, model::Structure& structure) {
    std::string line;
    std::shared_ptr<model::Residue> currentResidue = nullptr;
    bool success = true;
    
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        
        PdbParserStructures::RecordType recordType = PdbParserStructures::getRecordType(line);
        
        switch (recordType) {
            case PdbParserStructures::RecordType::ATOM:
            case PdbParserStructures::RecordType::HETATM:
                success = PdbParserRecords::parseAtomRecord(line, recordType, structure, currentResidue);
                break;
            case PdbParserStructures::RecordType::TER:
                success = PdbParserRecords::parseTerRecord(line, currentResidue, structure);
                break;
            case PdbParserStructures::RecordType::HELIX:
                success = PdbParserStructureRecords::parseHelixRecord(line, structure);
                break;
            case PdbParserStructures::RecordType::SHEET:
                success = PdbParserStructureRecords::parseSheetRecord(line, structure);
                break;
            case PdbParserStructures::RecordType::SSBOND:
                success = PdbParserStructureRecords::parseSSBondRecord(line, structure);
                break;
            case PdbParserStructures::RecordType::CRYST1:
                success = PdbParserStructureRecords::parseCryst1Record(line, structure);
                break;
            default:
                continue;
        }
        
        if (!success) {
            return false;
        }
    }
    
    // Calculate center of mass for the last residue if not already done
    if (currentResidue) {
        currentResidue->calculate_center_of_mass();
    }
    
    return true;
}

} // namespace structure
} // namespace io
} // namespace pygcmc