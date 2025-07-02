// src/io/forcefield/PrmParserDrudeScan.cpp

#include "PrmParserDrudeScan.hpp"
#include "PrmParserStructures.hpp"
#include "model/ModelModule.hpp"
#include <iostream>

namespace pygcmc {
namespace io {

void PrmParserDrudeScan::prescanForDrudeParameters(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output) {
    
    // First pass: collect all ATOM lines with ALPHA/THOLE from the entire file
    auto startPos = input.tellg();
    std::string scanLine;
    bool in_patch = false;  // Track if we're inside a PRES (patch) block
    int pres_count = 0;  // Count PRES blocks
    int atom_lines = 0;
    
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
}

} // namespace io
} // namespace pygcmc