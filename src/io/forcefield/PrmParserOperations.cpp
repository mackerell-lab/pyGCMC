// src/io/forcefield/PrmParserOperations.cpp

#include "PrmParserOperations.hpp"
#include "PrmParserSections.hpp"
#include "PrmParserStructures.hpp"
#include "PrmParserDrudeScan.hpp"
#include "PrmParserSectionHandlers.hpp"
#include "PrmParserBondedSections.hpp"
#include "PrmParserDrudeSections.hpp"
#include "model/ModelModule.hpp"
#include <fstream>

namespace pygcmc {
namespace io {

// Static debug flag shared across all parser components
static bool debug_output = false;

bool& PrmParserOperations::getDebugFlag() {
    return debug_output;
}

void PrmParserOperations::parseStream(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    bool inSection = false;
    std::string currentSection;
    
    // Sync debug flags
    PrmParserSections::getDebugFlag() = debug_output;
    
    // Debug output
    static int call_count = 0;
    call_count++;
    if (debug_output) std::cerr << "DEBUG: parseStream called #" << call_count << ", debug_output=" << debug_output << std::endl;
    int total_lines = 0;
    int atom_lines = 0;
    std::string first_line, last_line;
    
    // First pass: collect all ATOM lines with ALPHA/THOLE from the entire file
    PrmParserDrudeScan::prescanForDrudeParameters(input, ff, debug_output);
    
    // Now do normal parsing
    while (std::getline(input, line)) {
        total_lines++;
        if (total_lines == 1) first_line = line;
        last_line = line;
        if (debug_output) std::cerr << "Raw line: [" << line << "]" << std::endl;
        
        if (PrmParserStructures::isCommentLine(line)) {
            if (debug_output) std::cerr << "Skipping comment line" << std::endl;
            continue;
        }
        
        std::string cleanLine = PrmParserStructures::removeComments(line);
        cleanLine = PrmParserStructures::trim(cleanLine);
        if (debug_output) std::cerr << "Cleaned line: [" << cleanLine << "]" << std::endl;
        
        if (cleanLine.empty()) continue;
        
        // Skip topology lines from STR files (ATOM lines with ALPHA/THOLE were handled in pre-scan)
        if (PrmParserStructures::isTopologyLine(cleanLine)) {
            if (debug_output) std::cerr << "Skipping topology line: " << cleanLine << std::endl;
            continue;
        }
        
        // Note: In STR files, "end" may mark the end of topology section,
        // not the end of the file. Continue parsing for parameter sections.
        if (cleanLine == "END" || cleanLine == "end") {
            inSection = false;
            currentSection.clear();
            continue;
        }
        
        if (PrmParserStructures::isAtomsSection(cleanLine)) {
            if (debug_output) std::cerr << "Found ATOMS section" << std::endl;
            inSection = true;
            currentSection = "ATOMS";
            PrmParserSectionHandlers::parseAtomsSection(input, ff, debug_output);
        } else if (PrmParserStructures::isBondsSection(cleanLine)) {
            if (debug_output) std::cerr << "Found BONDS section" << std::endl;
            inSection = true;
            currentSection = "BONDS";
            PrmParserSectionHandlers::parseBondsSection(input, ff, debug_output);
        } else if (PrmParserStructures::isAnglesSection(cleanLine)) {
            if (debug_output) std::cerr << "Found ANGLES section" << std::endl;
            inSection = true;
            currentSection = "ANGLES";
            PrmParserSectionHandlers::parseAnglesSection(input, ff, debug_output);
        } else if (PrmParserStructures::isDihedralsSection(cleanLine)) {
            if (debug_output) std::cerr << "Found DIHEDRALS section" << std::endl;
            inSection = true;
            currentSection = "DIHEDRALS";
            PrmParserBondedSections::parseDihedralsSection(input, ff, debug_output);
        } else if (PrmParserStructures::isImproperSection(cleanLine)) {
            if (debug_output) std::cerr << "Found IMPROPER section" << std::endl;
            inSection = true;
            currentSection = "IMPROPER";
            PrmParserBondedSections::parseImproperSection(input, ff, debug_output);
        } else if (PrmParserStructures::isNonbondedSection(cleanLine)) {
            if (debug_output) std::cerr << "Found NONBONDED section" << std::endl;
            inSection = true;
            currentSection = "NONBONDED";
            
            // If the line starts with cutnb, we need to look for the previous NONBONDED line
            if (cleanLine.find("NONBONDED") == std::string::npos && cleanLine.find("cutnb") != std::string::npos) {
                // Store current position
                auto currentPos = input.tellg();
                std::string prevLine;
                
                // Go back to beginning of file
                input.seekg(0);
                
                // Read lines until we reach our current position
                while (input.tellg() < currentPos && std::getline(input, prevLine)) {
                    if (PrmParserStructures::isCommentLine(prevLine)) continue;
                    prevLine = PrmParserStructures::removeComments(prevLine);
                    prevLine = PrmParserStructures::trim(prevLine);
                    if (prevLine.find("NONBONDED") != std::string::npos) {
                        // Found the NONBONDED line, combine it with current line
                        cleanLine = prevLine + " " + cleanLine;
                        break;
                    }
                }
                
                // Restore position
                input.seekg(currentPos);
            }
            
            PrmParserSections::parseNonbondedSection(input, ff, cleanLine);
        } else if (PrmParserStructures::isNBFixSection(cleanLine)) {
            if (debug_output) std::cerr << "Found NBFIX section" << std::endl;
            inSection = true;
            currentSection = "NBFIX";
            PrmParserBondedSections::parseNBFixSection(input, ff, debug_output);
        } else if (PrmParserStructures::isAlphaTHoleSection(cleanLine)) {
            if (debug_output) std::cerr << "Found ALPHA/THOLE section" << std::endl;
            inSection = true;
            currentSection = "ALPHA";
            PrmParserDrudeSections::parseAlphaTHoleSection(input, ff, debug_output);
        } else if (PrmParserStructures::isLonePairSection(cleanLine)) {
            if (debug_output) std::cerr << "Found LONEPAIR section" << std::endl;
            inSection = true;
            currentSection = "LONEPAIR";
            PrmParserDrudeSections::parseLonePairSection(input, ff, debug_output);
        } else if (PrmParserStructures::isAnisotropySection(cleanLine)) {
            if (debug_output) std::cerr << "Found ANISOTROPY section" << std::endl;
            inSection = true;
            currentSection = "ANISOTROPY";
            PrmParserDrudeSections::parseAnisotropySection(input, ff, debug_output);
        } else if (!inSection) {
            // Handle non-section content if needed
            auto tokens = PrmParserStructures::tokenize(cleanLine);
            if (!tokens.empty()) {
                // Process any standalone parameters or commands
                if (tokens[0] == "set" || tokens[0] == "if" || tokens[0] == "read" || 
                    tokens[0] == "return" || tokens[0] == "BOMLEV" || tokens[0] == "WRNLEV") {
                    // Skip CHARMM control statements
                    continue;
                }
            }
        }
    }
    
    // Debug output
    if (debug_output) {
        std::cerr << "DEBUG: parseStream #" << call_count << " finished. Total lines: " << total_lines 
                  << ", ATOM lines with ALPHA/THOLE: " << atom_lines << std::endl;
        if (total_lines > 0) {
            std::cerr << "  First line: " << first_line << std::endl;
            std::cerr << "  Last line: " << last_line << std::endl;
        }
    }
}

} // namespace io
} // namespace pygcmc