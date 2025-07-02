// src/io/forcefield/PrmParserOperations.cpp

#include "PrmParserOperations.hpp"
#include "PrmParserSections.hpp"
#include "PrmParserStructures.hpp"
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
    auto startPos = input.tellg();
    std::string scanLine;
    bool in_patch = false;  // Track if we're inside a PRES (patch) block
    int pres_count = 0;  // Count PRES blocks
    
    while (std::getline(input, scanLine)) {
        if (PrmParserStructures::isCommentLine(scanLine)) continue;
        
        std::string cleanScanLine = PrmParserStructures::removeComments(scanLine);
        cleanScanLine = PrmParserStructures::trim(cleanScanLine);
        
        // Check if we're entering a PRES (patch) block - skip these
        if (cleanScanLine.find("PRES ") == 0) {
            in_patch = true;
            pres_count++;
            if (debug_output) {
                std::cerr << "Pre-scan: Entering PRES block #" << pres_count << ": " << cleanScanLine << std::endl;
            }
            continue;
        }
        
        
        // Check if we're exiting a patch block
        if (in_patch && (cleanScanLine == "end" || cleanScanLine == "END")) {
            in_patch = false;
            if (debug_output) {
                std::cerr << "Pre-scan: Exiting PRES block" << std::endl;
            }
            continue;
        }
        
        // For ATOM lines with ALPHA/THOLE, skip if inside patch blocks
        // But LONEPAIR and ANISOTROPY should be processed even in patch blocks
        
        if (cleanScanLine.find("ATOM ") == 0 &&
            (cleanScanLine.find("ALPHA") != std::string::npos || cleanScanLine.find("THOLE") != std::string::npos))
        {
            // Skip ATOM lines with ALPHA/THOLE if inside patch blocks
            if (in_patch) {
                if (debug_output) {
                    std::cerr << "Pre-scan: Skipping ATOM line in PRES block: " << cleanScanLine << std::endl;
                }
                continue;
            }
            // Parse ALPHA/THOLE from ATOM line
            auto tokens = PrmParserStructures::tokenize(cleanScanLine);
            if (tokens.size() >= 4) {
                std::string atomType = tokens[2];  // The atom type (e.g., ODW)
                
                // Note: For CHARMM force fields, we keep the last occurrence, not the first
                // So we don't skip duplicates here
                
                // Special logging for ND2A2
                if (atomType == "ND2A2" && debug_output) {
                    std::cerr << "Pre-scan: Processing ND2A2 line: " << cleanScanLine << std::endl;
                    std::cerr << "Pre-scan: in_patch = " << in_patch << std::endl;
                }
                
                // Find ALPHA and THOLE values
                for (size_t i = 0; i < tokens.size(); ++i) {
                    if (tokens[i] == "ALPHA" && i + 1 < tokens.size()) {
                        try {
                            double alpha = PrmParserStructures::safe_stod(tokens[i + 1], "ALPHA for " + atomType);
                            double thole = 0.0; // Default thole
                            
                            // Look for THOLE after ALPHA
                            if (i + 3 < tokens.size() && tokens[i + 2] == "THOLE") {
                                thole = PrmParserStructures::safe_stod(tokens[i + 3], "THOLE for " + atomType);
                            }
                            
                            ff.add_alpha_thole_params(atomType, alpha, thole);
                            atom_lines++;
                            if (debug_output) {
                                std::cerr << "Pre-scan: Added ALPHA/THOLE for " << atomType
                                          << ": alpha=" << alpha << ", thole=" << thole << std::endl;
                            }
                            break; 
                        } catch (const std::exception& e) {
                            if (debug_output) {
                                std::cerr << "Pre-scan: Failed to parse ALPHA/THOLE from line: " << cleanScanLine << std::endl;
                                std::cerr << "Error: " << e.what() << std::endl;
                            }
                        }
                    }
                }
            }
        }
        
        // Check for LONEPAIR lines (include those in PRES blocks)
        if (cleanScanLine.find("LONEPAIR") == 0) {
            auto tokens = PrmParserStructures::tokenize(cleanScanLine);
            if (tokens.size() >= 8) {
                try {
                    pygcmc::model::LonePairParams params;
                    params.type = tokens[1];  // e.g., "bisector" or "relative"
                    params.atom1 = tokens[2]; // The lonepair atom
                    params.host = tokens[3];  // Host atom
                    params.atom2 = tokens[4]; // Reference atom 1
                    params.atom3 = tokens[5]; // Reference atom 2
                    
                    // Find "distance", "angle", "dihe" keywords
                    for (size_t i = 6; i < tokens.size(); ++i) {
                        if (tokens[i] == "distance" && i + 1 < tokens.size()) {
                            params.distance = PrmParserStructures::safe_stod(tokens[i + 1], "LONEPAIR distance");
                        } else if (tokens[i] == "angle" && i + 1 < tokens.size()) {
                            params.angle = PrmParserStructures::safe_stod(tokens[i + 1], "LONEPAIR angle");
                        } else if (tokens[i] == "dihe" && i + 1 < tokens.size()) {
                            params.dihedral = PrmParserStructures::safe_stod(tokens[i + 1], "LONEPAIR dihedral");
                        }
                    }
                    
                    ff.add_lonepair(params);
                    if (debug_output) {
                        std::cerr << "Pre-scan: Added LONEPAIR " << params.type << " for " << params.atom1 << std::endl;
                    }
                } catch (const std::exception& e) {
                    if (debug_output) {
                        std::cerr << "Pre-scan: Failed to parse LONEPAIR from line: " << cleanScanLine << std::endl;
                        std::cerr << "Error: " << e.what() << std::endl;
                    }
                }
            }
        }
        
        // Check for ANISOTROPY lines (include those in PRES blocks)
        if (cleanScanLine.find("ANISOTROPY") == 0) {
            auto tokens = PrmParserStructures::tokenize(cleanScanLine);
            if (tokens.size() >= 8) {
                try {
                    pygcmc::model::AnisotropyParams params;
                    params.type = tokens[1];  // atom type
                    // tokens[2], [3], [4] are reference atoms, but we'll skip them for now
                    
                    // Find A11, A22, A33 values
                    for (size_t i = 5; i < tokens.size(); ++i) {
                        if (tokens[i] == "A11" && i + 1 < tokens.size()) {
                            params.a11 = PrmParserStructures::safe_stod(tokens[i + 1], "ANISOTROPY A11");
                        } else if (tokens[i] == "A22" && i + 1 < tokens.size()) {
                            params.a22 = PrmParserStructures::safe_stod(tokens[i + 1], "ANISOTROPY A22");
                        } else if (tokens[i] == "A33" && i + 1 < tokens.size()) {
                            params.a33 = PrmParserStructures::safe_stod(tokens[i + 1], "ANISOTROPY A33");
                        }
                    }
                    
                    ff.add_anisotropy(params);
                    if (debug_output) {
                        std::cerr << "Pre-scan: Added ANISOTROPY for " << params.type << std::endl;
                    }
                } catch (const std::exception& e) {
                    if (debug_output) {
                        std::cerr << "Pre-scan: Failed to parse ANISOTROPY from line: " << cleanScanLine << std::endl;
                        std::cerr << "Error: " << e.what() << std::endl;
                    }
                }
            }
        }
        
        // Check for THOLE global parameters line (THOLE TCUT ... MAXNBTHOLE ...)
        if (cleanScanLine.find("THOLE") == 0 && cleanScanLine.find("TCUT") != std::string::npos) {
            auto tokens = PrmParserStructures::tokenize(cleanScanLine);
            double tcut = 5.0;  // default
            int maxnbthole = 5000;  // default
            
            for (size_t i = 0; i < tokens.size(); ++i) {
                if (tokens[i] == "TCUT" && i + 1 < tokens.size()) {
                    tcut = PrmParserStructures::safe_stod(tokens[i + 1], "TCUT");
                } else if (tokens[i] == "MAXNBTHOLE" && i + 1 < tokens.size()) {
                    maxnbthole = PrmParserStructures::safe_stoi(tokens[i + 1], "MAXNBTHOLE");
                }
            }
            
            ff.set_drude_global_params(tcut, maxnbthole);
            if (debug_output) {
                std::cerr << "Pre-scan: Set THOLE global params: TCUT=" << tcut 
                          << ", MAXNBTHOLE=" << maxnbthole << std::endl;
            }
        }
        
        // Check for NBTHOLE lines (atom1 atom2 value)
        // These appear after the THOLE global line and before "end"
        if (!cleanScanLine.empty() && !in_patch) {
            auto tokens = PrmParserStructures::tokenize(cleanScanLine);
            if (tokens.size() == 3) {
                // Check if this might be an NBTHOLE line (two atom types + value)
                // NBTHOLE lines have format: ATOM1 ATOM2 VALUE
                // Atom types typically contain letters and numbers
                bool isAtomType1 = tokens[0].find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") == std::string::npos;
                bool isAtomType2 = tokens[1].find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789") == std::string::npos;
                
                if (isAtomType1 && isAtomType2 && tokens[0].length() <= 8 && tokens[1].length() <= 8) {
                    try {
                        // Try to parse the third token as a double
                        double thole = PrmParserStructures::safe_stod(tokens[2], "potential NBTHOLE value");
                        
                        // Additional checks to avoid parsing other sections
                        if (tokens[0] != "MASS" && tokens[0] != "ATOM" && tokens[0] != "BOND" &&
                            tokens[0] != "ANGLE" && tokens[0] != "DIHE" && tokens[0] != "IMPR" &&
                            tokens[0] != "NONB" && tokens[0] != "NBFIX" && tokens[0] != "HBOND" &&
                            tokens[0] != "CMAP" && tokens[0] != "END" && tokens[0] != "end") {
                            // Likely an NBTHOLE line
                            ff.add_nbthole(tokens[0], tokens[1], thole);
                            if (debug_output) {
                                std::cerr << "Pre-scan: Added NBTHOLE for " << tokens[0] 
                                          << "-" << tokens[1] << " = " << thole << std::endl;
                            }
                        }
                    } catch (...) {
                        // Not an NBTHOLE line, ignore
                    }
                }
            }
        }
    }
    
    if (debug_output) {
        std::cerr << "Pre-scan complete: found " << atom_lines << " ATOM lines with ALPHA/THOLE" << std::endl;
        std::cerr << "Pre-scan: Detected " << pres_count << " PRES blocks" << std::endl;
        if (in_patch) {
            std::cerr << "WARNING: Pre-scan ended while still in_patch=true" << std::endl;
        }
    }
    
    // Reset to beginning of stream for normal parsing
    input.clear();
    input.seekg(startPos);
    
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
            parseAtomsSection(input, ff);
        } else if (PrmParserStructures::isBondsSection(cleanLine)) {
            if (debug_output) std::cerr << "Found BONDS section" << std::endl;
            inSection = true;
            currentSection = "BONDS";
            parseBondsSection(input, ff);
        } else if (PrmParserStructures::isAnglesSection(cleanLine)) {
            if (debug_output) std::cerr << "Found ANGLES section" << std::endl;
            inSection = true;
            currentSection = "ANGLES";
            parseAnglesSection(input, ff);
        } else if (PrmParserStructures::isDihedralsSection(cleanLine)) {
            if (debug_output) std::cerr << "Found DIHEDRALS section" << std::endl;
            inSection = true;
            currentSection = "DIHEDRALS";
            PrmParserSections::parseDihedralsSection(input, ff);
        } else if (PrmParserStructures::isImproperSection(cleanLine)) {
            if (debug_output) std::cerr << "Found IMPROPER section" << std::endl;
            inSection = true;
            currentSection = "IMPROPER";
            PrmParserSections::parseImproperSection(input, ff);
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
            PrmParserSections::parseNBFixSection(input, ff);
        } else if (PrmParserStructures::isAlphaTHoleSection(cleanLine)) {
            if (debug_output) std::cerr << "Found ALPHA/THOLE section" << std::endl;
            inSection = true;
            currentSection = "ALPHA";
            PrmParserSections::parseAlphaTHoleSection(input, ff);
        } else if (PrmParserStructures::isLonePairSection(cleanLine)) {
            if (debug_output) std::cerr << "Found LONEPAIR section" << std::endl;
            inSection = true;
            currentSection = "LONEPAIR";
            PrmParserSections::parseLonePairSection(input, ff);
        } else if (PrmParserStructures::isAnisotropySection(cleanLine)) {
            if (debug_output) std::cerr << "Found ANISOTROPY section" << std::endl;
            inSection = true;
            currentSection = "ANISOTROPY";
            PrmParserSections::parseAnisotropySection(input, ff);
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

void PrmParserOperations::parseAtomsSection(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    if (debug_output) std::cerr << "\n=== Entering ATOMS/MASS section parsing ===" << std::endl;
    
    // Add debug counter
    int lines_processed = 0;
    int atoms_with_alpha = 0;
    
    while (std::getline(input, line)) {
        lines_processed++;
        if (debug_output) std::cerr << "Raw atom line: [" << line << "]" << std::endl;
        
        if (PrmParserStructures::isCommentLine(line)) {
            if (debug_output) std::cerr << "Skipping comment line in atoms section" << std::endl;
            continue;
        }
        
        std::string fullLine = PrmParserStructures::readContinuationLine(input, line);
        fullLine = PrmParserStructures::removeComments(fullLine);
        fullLine = PrmParserStructures::trim(fullLine);
        
        if (fullLine.empty()) {
            if (debug_output) std::cerr << "Skipping empty line in atoms section" << std::endl;
            continue;
        }
        
        if (debug_output) std::cerr << "Processed line: [" << fullLine << "]" << std::endl;
        
        // Note: ATOM lines with ALPHA/THOLE are handled in the pre-scan phase
        // parseAtomsSection should only parse MASS entries, not ATOM lines
        
        // Skip other topology lines from STR files
        if (PrmParserStructures::isTopologyLine(fullLine)) {
            if (debug_output) std::cerr << "Skipping topology line in atoms section: " << fullLine << std::endl;
            continue;
        }
        
        // Check for section end
        if (fullLine == "END" || fullLine == "end" || PrmParserStructures::isBondsSection(fullLine) || PrmParserStructures::isAnglesSection(fullLine) || 
            PrmParserStructures::isDihedralsSection(fullLine) || PrmParserStructures::isImproperSection(fullLine) || 
            PrmParserStructures::isNonbondedSection(fullLine) || PrmParserStructures::isNBFixSection(fullLine)) {
            if (debug_output) std::cerr << "Found section end marker: " << fullLine << std::endl;
            // Note: For atoms section, we use the original line length, not fullLine
            // because fullLine may be concatenated from multiple lines
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(fullLine);
        if (debug_output) {
            std::cerr << "Tokens:";
            for (const auto& token : tokens) {
                std::cerr << " [" << token << "]";
            }
            std::cerr << std::endl;
        }
        
        if (tokens.size() >= 4 && tokens[0] == "MASS") {
            std::string atomType = tokens[2];
            double mass = PrmParserStructures::safe_stod(tokens[3], "atom mass for type " + atomType);
            ff.add_atom_mass(atomType, mass);
            if (debug_output) std::cerr << "Added atom mass: " << atomType << " = " << mass << std::endl;
        } else {
            if (debug_output) std::cerr << "Skipping line: not a valid MASS entry" << std::endl;
        }
    }
    
    if (debug_output) {
        std::cerr << "\n=== Finished ATOMS/MASS section parsing ===" << std::endl;
        std::cerr << "Final atom_masses size: " << ff.get_num_atom_types() << std::endl;
        std::cerr << "\nStored atom masses:" << std::endl;
        for (const auto& pair : ff.get_atom_masses()) {
            std::cerr << pair.first << " = " << pair.second << std::endl;
        }
    }
    
    // Add debug summary
    if (debug_output) {
        std::cerr << "\n=== Exiting ATOMS/MASS section parsing ===" << std::endl;
        std::cerr << "Lines processed: " << lines_processed << std::endl;
        std::cerr << "ATOM lines with ALPHA/THOLE found: " << atoms_with_alpha << std::endl;
    }
}

void PrmParserOperations::parseBondsSection(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    if (debug_output) std::cerr << "\n=== Entering BONDS section parsing ===" << std::endl;
    
    while (std::getline(input, line)) {
        if (debug_output) std::cerr << "Raw bond line: [" << line << "]" << std::endl;
        
        // Save original line length before any modifications
        size_t originalLineLength = line.length();
        
        if (PrmParserStructures::isCommentLine(line)) {
            if (debug_output) std::cerr << "Skipping comment line in bonds section" << std::endl;
            continue;
        }
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) {
            if (debug_output) std::cerr << "Skipping empty line in bonds section" << std::endl;
            continue;
        }
        
        // Skip topology lines from STR files
        if (PrmParserStructures::isTopologyLine(line)) {
            if (debug_output) std::cerr << "Skipping topology line in bonds section: " << line << std::endl;
            continue;
        }
        
        // Check for section end
        if (line == "END" || line == "end" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isAnglesSection(line) || 
            PrmParserStructures::isDihedralsSection(line) || PrmParserStructures::isImproperSection(line) || 
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line)) {
            if (debug_output) std::cerr << "Found section end marker in bonds: " << line << std::endl;
            input.seekg(-static_cast<std::streamoff>(originalLineLength + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        if (debug_output) {
            std::cerr << "Bond tokens:";
            for (const auto& token : tokens) {
                std::cerr << " [" << token << "]";
            }
            std::cerr << std::endl;
        }
        
        if (tokens.size() >= 4) {
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            double kb = PrmParserStructures::safe_stod(tokens[2], "bond Kb for " + type1 + "-" + type2);
            double b0 = PrmParserStructures::safe_stod(tokens[3], "bond b0 for " + type1 + "-" + type2);
            
            if (debug_output) {
                std::cerr << "\nProcessing bond: " << type1 << "-" << type2 << std::endl;
                std::cerr << "  kb = " << kb << ", b0 = " << b0 << std::endl;
            }
            
            ff.add_bond_params(type1, type2, kb, b0);
            
            if (debug_output) {
                std::cerr << "  Current bond_params size: " << ff.get_num_bond_types() << std::endl;
            }
        }
    }
    
    if (debug_output) {
        std::cerr << "\n=== Finished BONDS section parsing ===" << std::endl;
        std::cerr << "Final bond_params size: " << ff.get_num_bond_types() << std::endl;
    }
}

void PrmParserOperations::parseAnglesSection(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    if (debug_output) std::cerr << "\n=== Entering ANGLES section parsing ===" << std::endl;
    
    while (std::getline(input, line)) {
        if (debug_output && line.find("ATOM") != std::string::npos) {
            std::cerr << "WARNING: ATOM line in angles section: [" << line << "]" << std::endl;
        }
        if (debug_output) std::cerr << "Raw angle line: [" << line << "]" << std::endl;
        
        // Save original line length before any modifications
        size_t originalLineLength = line.length();
        
        if (PrmParserStructures::isCommentLine(line)) {
            if (debug_output) std::cerr << "Skipping comment line in angles section" << std::endl;
            continue;
        }
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) continue;
        
        // Skip ATOM lines (they're not angle parameters)
        if (line.find("ATOM ") == 0) {
            if (debug_output) std::cerr << "Skipping ATOM line in angles section: " << line << std::endl;
            continue;
        }
        
        // Skip topology lines from STR files
        if (PrmParserStructures::isTopologyLine(line)) {
            if (debug_output) std::cerr << "Skipping topology line in angles section: " << line << std::endl;
            continue;
        }
        
        // Check for section end
        if (line == "END" || line == "end" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isBondsSection(line) || 
            PrmParserStructures::isDihedralsSection(line) || PrmParserStructures::isImproperSection(line) || 
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line)) {
            input.seekg(-static_cast<std::streamoff>(originalLineLength + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        if (tokens.size() >= 5) {
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            std::string type3 = tokens[2];
            double ktheta = PrmParserStructures::safe_stod(tokens[3], "angle Ktheta for " + type1 + "-" + type2 + "-" + type3);
            double theta0 = PrmParserStructures::safe_stod(tokens[4], "angle theta0 for " + type1 + "-" + type2 + "-" + type3);
            
            double kub = 0.0;
            double s0 = 0.0;
            if (tokens.size() >= 7) {
                kub = PrmParserStructures::safe_stod(tokens[5], "angle Kub for " + type1 + "-" + type2 + "-" + type3);
                s0 = PrmParserStructures::safe_stod(tokens[6], "angle S0 for " + type1 + "-" + type2 + "-" + type3);
            }
            
            if (debug_output) {
                std::cerr << "\nProcessing angle: " << type1 << "-" << type2 << "-" << type3 << std::endl;
                std::cerr << "  ktheta = " << ktheta << ", theta0 = " << theta0 << std::endl;
            }
            
            ff.add_angle_params(type1, type2, type3, ktheta, theta0, kub, s0);
            
            if (debug_output) {
                std::cerr << "  Current angle_params size: " << ff.get_num_angle_types() << std::endl;
            }
        }
    }
    
    if (debug_output) {
        std::cerr << "\n=== Finished ANGLES section parsing ===" << std::endl;
        std::cerr << "Final angle_params size: " << ff.get_num_angle_types() << std::endl;
    }
}

} // namespace io
} // namespace pygcmc