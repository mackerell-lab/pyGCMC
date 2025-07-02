// src/io/forcefield/PrmParserDrudeSections.cpp

#include "PrmParserDrudeSections.hpp"
#include "PrmParserStructures.hpp"
#include "model/ModelModule.hpp"

namespace pygcmc {
namespace io {

void PrmParserDrudeSections::parseAlphaTHoleSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output) {
    std::string line;
    while (std::getline(input, line)) {
        // Save original line length before any modifications
        size_t originalLineLength = line.length();
        
        if (PrmParserStructures::isCommentLine(line)) continue;
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) continue;
        
        // Skip topology lines from STR files
        if (PrmParserStructures::isTopologyLine(line)) {
            if (debug_output) std::cerr << "Skipping topology line in ALPHA/THOLE section: " << line << std::endl;
            continue;
        }
        
        // Check for section end
        if (line == "END" || line == "end" || PrmParserStructures::isAtomsSection(line) || 
            PrmParserStructures::isBondsSection(line) || PrmParserStructures::isAnglesSection(line) || 
            PrmParserStructures::isDihedralsSection(line) || PrmParserStructures::isImproperSection(line) ||
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line) ||
            PrmParserStructures::isLonePairSection(line) || PrmParserStructures::isAnisotropySection(line)) {
            input.seekg(-static_cast<std::streamoff>(originalLineLength + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        // Format: atomType  alpha  thole
        if (tokens.size() >= 3) {
            try {
                std::string type = tokens[0];
                double alpha = PrmParserStructures::safe_stod(tokens[1], "alpha for " + type);
                double thole = PrmParserStructures::safe_stod(tokens[2], "thole for " + type);
                
                ff.add_alpha_thole_params(type, alpha, thole);
                
                if (debug_output) {
                    std::cerr << "Added ALPHA/THOLE for " << type 
                              << ": alpha = " << alpha << ", thole = " << thole << std::endl;
                }
            } catch (const std::exception& e) {
                if (debug_output) std::cerr << "Warning: Skipping ALPHA/THOLE line due to parsing error: " << line << std::endl;
                continue;
            }
        }
    }
}

void PrmParserDrudeSections::parseLonePairSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output) {
    std::string line;
    while (std::getline(input, line)) {
        // Save original line length before any modifications
        size_t originalLineLength = line.length();
        
        if (PrmParserStructures::isCommentLine(line)) continue;
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) continue;
        
        // Skip topology lines from STR files
        if (PrmParserStructures::isTopologyLine(line)) {
            if (debug_output) std::cerr << "Skipping topology line in LONEPAIR section: " << line << std::endl;
            continue;
        }
        
        // Check for section end
        if (line == "END" || line == "end" || PrmParserStructures::isAtomsSection(line) || 
            PrmParserStructures::isBondsSection(line) || PrmParserStructures::isAnglesSection(line) || 
            PrmParserStructures::isDihedralsSection(line) || PrmParserStructures::isImproperSection(line) ||
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line) ||
            PrmParserStructures::isAlphaTHoleSection(line) || PrmParserStructures::isAnisotropySection(line)) {
            input.seekg(-static_cast<std::streamoff>(originalLineLength + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        // Various formats for LONEPAIR definitions
        // Common format: type host atom1 atom2 atom3 distance angle dihedral
        if (tokens.size() >= 5) {
            try {
                pygcmc::model::LonePairParams params;
                params.type = tokens[0];
                params.host = tokens[1];
                params.atom1 = tokens[2];
                params.atom2 = tokens[3];
                
                // Handle different lonepair formats
                if (tokens.size() >= 6) {
                    params.atom3 = tokens[4];
                    if (tokens.size() >= 7) {
                        params.distance = PrmParserStructures::safe_stod(tokens[5], "lonepair distance");
                        if (tokens.size() >= 8) {
                            params.angle = PrmParserStructures::safe_stod(tokens[6], "lonepair angle");
                            if (tokens.size() >= 9) {
                                params.dihedral = PrmParserStructures::safe_stod(tokens[7], "lonepair dihedral");
                            }
                        }
                    }
                } else {
                    // Some formats have only 4 atoms
                    params.atom3 = tokens[4];
                }
                
                ff.add_lonepair(params);
                
                if (debug_output) {
                    std::cerr << "Added LONEPAIR: type=" << params.type 
                              << ", host=" << params.host << std::endl;
                }
            } catch (const std::exception& e) {
                if (debug_output) std::cerr << "Warning: Skipping LONEPAIR line due to parsing error: " << line << std::endl;
                continue;
            }
        }
    }
}

void PrmParserDrudeSections::parseAnisotropySection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output) {
    std::string line;
    while (std::getline(input, line)) {
        // Save original line length before any modifications
        size_t originalLineLength = line.length();
        
        if (PrmParserStructures::isCommentLine(line)) continue;
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) continue;
        
        // Skip topology lines from STR files
        if (PrmParserStructures::isTopologyLine(line)) {
            if (debug_output) std::cerr << "Skipping topology line in ANISOTROPY section: " << line << std::endl;
            continue;
        }
        
        // Check for section end
        if (line == "END" || line == "end" || PrmParserStructures::isAtomsSection(line) || 
            PrmParserStructures::isBondsSection(line) || PrmParserStructures::isAnglesSection(line) || 
            PrmParserStructures::isDihedralsSection(line) || PrmParserStructures::isImproperSection(line) ||
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line) ||
            PrmParserStructures::isAlphaTHoleSection(line) || PrmParserStructures::isLonePairSection(line)) {
            input.seekg(-static_cast<std::streamoff>(originalLineLength + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        // Format: atomType a11 a22 [a33]
        if (tokens.size() >= 3) {
            try {
                pygcmc::model::AnisotropyParams params;
                params.type = tokens[0];
                params.a11 = PrmParserStructures::safe_stod(tokens[1], "anisotropy a11 for " + params.type);
                params.a22 = PrmParserStructures::safe_stod(tokens[2], "anisotropy a22 for " + params.type);
                
                // Optional a33 component
                if (tokens.size() >= 4) {
                    params.a33 = PrmParserStructures::safe_stod(tokens[3], "anisotropy a33 for " + params.type);
                }
                
                ff.add_anisotropy(params);
                
                if (debug_output) {
                    std::cerr << "Added ANISOTROPY for " << params.type 
                              << ": a11=" << params.a11 << ", a22=" << params.a22;
                    if (tokens.size() >= 4) {
                        std::cerr << ", a33=" << params.a33;
                    }
                    std::cerr << std::endl;
                }
            } catch (const std::exception& e) {
                if (debug_output) std::cerr << "Warning: Skipping ANISOTROPY line due to parsing error: " << line << std::endl;
                continue;
            }
        }
    }
}

} // namespace io
} // namespace pygcmc