// src/io/forcefield/PrmParserSections.cpp

#include "PrmParserSections.hpp"
#include "PrmParserStructures.hpp"
#include "PrmParserBondedSections.hpp"
#include "model/ModelModule.hpp"

namespace pygcmc {
namespace io {

// Static debug flag shared across all parser components
static bool debug_output = false;

bool& PrmParserSections::getDebugFlag() {
    return debug_output;
}

void PrmParserSections::parseNonbondedSection(std::istream& input, pygcmc::model::ForceField& ff, const std::string& firstLine) {
    std::string line = firstLine;
    std::string fullLine = PrmParserStructures::readContinuationLine(input, line);
    
    if (debug_output) std::cerr << "=== Entering NONBONDED section parsing ===" << std::endl;
    if (debug_output) std::cerr << "First line: [" << firstLine << "]" << std::endl;
    
    // Parse header parameters
    auto tokens = PrmParserStructures::tokenize(fullLine);
    if (!tokens.empty()) {
        if (debug_output) {
            std::cerr << "Initial tokens:";
            for (const auto& token : tokens) {
                std::cerr << " [" << token << "]";
            }
            std::cerr << std::endl;
        }
        
        // Skip the NONBONDED keyword
        if (tokens[0] == "NONBONDED") {
            if (debug_output) std::cerr << "Skipping NONBONDED keyword" << std::endl;
            tokens.erase(tokens.begin());
        }
        
        // Process parameters
        auto& params = ff.get_nonbonded_params();
        for (size_t i = 0; i < tokens.size(); ++i) {
            if (debug_output) std::cerr << "Processing parameter token: [" << tokens[i] << "]" << std::endl;
            
            if (tokens[i] == "nbxmod" && i + 1 < tokens.size()) {
                params.nbxmod = PrmParserStructures::safe_stoi(tokens[++i], "nbxmod");
                if (debug_output) std::cerr << "Set nbxmod = " << params.nbxmod << std::endl;
            } else if (tokens[i] == "cutnb" && i + 1 < tokens.size()) {
                params.cutnb = PrmParserStructures::safe_stod(tokens[++i], "cutnb");
                if (debug_output) std::cerr << "Set cutnb = " << params.cutnb << std::endl;
            } else if (tokens[i] == "ctofnb" && i + 1 < tokens.size()) {
                params.ctofnb = PrmParserStructures::safe_stod(tokens[++i], "ctofnb");
                if (debug_output) std::cerr << "Set ctofnb = " << params.ctofnb << std::endl;
            } else if (tokens[i] == "ctonnb" && i + 1 < tokens.size()) {
                params.ctonnb = PrmParserStructures::safe_stod(tokens[++i], "ctonnb");
                if (debug_output) std::cerr << "Set ctonnb = " << params.ctonnb << std::endl;
            } else if (tokens[i] == "eps" && i + 1 < tokens.size()) {
                params.eps = PrmParserStructures::safe_stod(tokens[++i], "eps");
                if (debug_output) std::cerr << "Set eps = " << params.eps << std::endl;
            } else if (tokens[i] == "e14fac" && i + 1 < tokens.size()) {
                params.e14fac = PrmParserStructures::safe_stod(tokens[++i], "e14fac");
                if (debug_output) std::cerr << "Set e14fac = " << params.e14fac << std::endl;
            } else if (tokens[i] == "wmin" && i + 1 < tokens.size()) {
                params.wmin = PrmParserStructures::safe_stod(tokens[++i], "wmin");
                if (debug_output) std::cerr << "Set wmin = " << params.wmin << std::endl;
            } else if (tokens[i] == "cdiel") {
                params.cdiel = true;
                if (debug_output) std::cerr << "Set cdiel = true" << std::endl;
            } else if (tokens[i] == "fshift") {
                params.fshift = true;
                if (debug_output) std::cerr << "Set fshift = true" << std::endl;
            } else if (tokens[i] == "vatom") {
                params.vatom = true;
                if (debug_output) std::cerr << "Set vatom = true" << std::endl;
            } else if (tokens[i] == "vdistance") {
                params.vdistance = true;
                if (debug_output) std::cerr << "Set vdistance = true" << std::endl;
            } else if (tokens[i] == "vfswitch") {
                params.vfswitch = true;
                if (debug_output) std::cerr << "Set vfswitch = true" << std::endl;
            }
        }
    }
    
    if (debug_output) std::cerr << "\n=== Starting atom type parameters parsing ===" << std::endl;
    if (debug_output) std::cerr << "Current lj_params map size: " << ff.get_num_lj_params() << std::endl;
    
    // Parse atom type parameters
    while (true) {
        const std::streampos linePos = input.tellg();
        if (!std::getline(input, line)) {
            break;
        }
        if (PrmParserStructures::isCommentLine(line)) continue;
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) {
            if (debug_output) std::cerr << "Skipping empty line" << std::endl;
            continue;
        }
        
        // Skip topology lines from STR files
        if (PrmParserStructures::isTopologyLine(line)) {
            if (debug_output) std::cerr << "Skipping topology line in nonbonded section: " << line << std::endl;
            continue;
        }
        
        if (debug_output) std::cerr << "Processing cleaned line: [" << line << "]" << std::endl;
        tokens = PrmParserStructures::tokenize(line);
        if (debug_output) {
            std::cerr << "Tokens:";
            for (const auto& token : tokens) {
                std::cerr << " [" << token << "]";
            }
            std::cerr << std::endl;
        }
        
        if (tokens.empty()) continue;
        
        // Check for section end or new section
        if (line == "END" || line == "end" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isBondsSection(line) || 
            PrmParserStructures::isAnglesSection(line) || PrmParserStructures::isDihedralsSection(line) || 
            PrmParserStructures::isImproperSection(line)) {
            if (debug_output) std::cerr << "Found section end marker: " << tokens[0] << std::endl;
            break;
        }

        // Drude-related sections can follow NONBONDED in CHARMM parameter files (e.g. "ALPHA").
        // Stop parsing NONBONDED and let the outer parser handle the new section header.
        if (PrmParserStructures::isAlphaTHoleSection(line) ||
            PrmParserStructures::isLonePairSection(line) ||
            PrmParserStructures::isAnisotropySection(line)) {
            input.clear();
            input.seekg(linePos);
            break;
        }
        
        // If we find NBFIX, parse it as a new section
        if (PrmParserStructures::isNBFixSection(line)) {
            if (debug_output) std::cerr << "Found NBFIX section" << std::endl;
            PrmParserBondedSections::parseNBFixSection(input, ff, debug_output);
            break;
        }
        
        // Format: atomType ignored epsilon Rmin [ignored ignored ignored ignored]
        if (tokens.size() >= 4) {
            std::string atomType = tokens[0];
            // Skip the ignored value (usually 0.0)
            try {
                double epsilon = PrmParserStructures::safe_stod(tokens[2], "LJ epsilon for " + atomType);
                double rmin_half = PrmParserStructures::safe_stod(tokens[3], "LJ Rmin/2 for " + atomType);
                
                if (debug_output) {
                    std::cerr << "\n*** Parsing atom type: " << atomType << " ***" << std::endl;
                    std::cerr << "  epsilon = " << epsilon << std::endl;
                    std::cerr << "  rmin_half = " << rmin_half << std::endl;
                }
                
                ff.add_lj_params(atomType, epsilon, rmin_half);
                
                if (debug_output) {
                    std::cerr << "Successfully stored " << atomType << " parameters" << std::endl;
                    std::cerr << "Current lj_params map size: " << ff.get_num_lj_params() << std::endl;
                }
            } catch (const std::exception& e) {
                throw std::runtime_error("Failed to parse LJ parameters for " + atomType + ": " + e.what());
            }
        } else if (!tokens.empty() && !PrmParserStructures::isCommentLine(line)) {
            throw std::runtime_error("Malformed NONBONDED parameters in line: " + line + 
                                   "\nExpected at least 4 tokens, got " + std::to_string(tokens.size()));
        }
    }
    
    if (debug_output) std::cerr << "\n=== Finished NONBONDED section parsing ===" << std::endl;
    if (debug_output) std::cerr << "Final lj_params map size: " << ff.get_num_lj_params() << std::endl;
}

} // namespace io
} // namespace pygcmc
