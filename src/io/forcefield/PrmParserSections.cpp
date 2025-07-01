// src/io/forcefield/PrmParserSections.cpp

#include "PrmParserSections.hpp"
#include "PrmParserStructures.hpp"
#include "model/ModelModule.hpp"

namespace pygcmc {
namespace io {

// Static debug flag shared across all parser components
static bool debug_output = false;

bool& PrmParserSections::getDebugFlag() {
    return debug_output;
}

void PrmParserSections::parseDihedralsSection(std::istream& input, pygcmc::model::ForceField& ff) {
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
            if (debug_output) std::cerr << "Skipping topology line in dihedrals section: " << line << std::endl;
            continue;
        }
        
        // Check for section end
        if (line == "END" || line == "end" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isBondsSection(line) || 
            PrmParserStructures::isAnglesSection(line) || PrmParserStructures::isImproperSection(line) || 
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line)) {
            input.seekg(-static_cast<std::streamoff>(originalLineLength + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        if (tokens.size() >= 7) {
            std::string type1 = tokens[0];
            std::string type2 = tokens[1];
            std::string type3 = tokens[2];
            std::string type4 = tokens[3];
            double kchi = PrmParserStructures::safe_stod(tokens[4], "dihedral Kchi for " + type1 + "-" + type2 + "-" + type3 + "-" + type4);
            int n = PrmParserStructures::safe_stoi(tokens[5], "dihedral n for " + type1 + "-" + type2 + "-" + type3 + "-" + type4);
            double delta = PrmParserStructures::safe_stod(tokens[6], "dihedral delta for " + type1 + "-" + type2 + "-" + type3 + "-" + type4);
            
            ff.add_dihedral_params(type1, type2, type3, type4, kchi, n, delta);
        }
    }
}

void PrmParserSections::parseImproperSection(std::istream& input, pygcmc::model::ForceField& ff) {
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
            if (debug_output) std::cerr << "Skipping topology line in improper section: " << line << std::endl;
            continue;
        }
        
        // Check for section end
        if (line == "END" || line == "end" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isBondsSection(line) || 
            PrmParserStructures::isAnglesSection(line) || PrmParserStructures::isDihedralsSection(line) || 
            PrmParserStructures::isNonbondedSection(line) || PrmParserStructures::isNBFixSection(line)) {
            input.seekg(-static_cast<std::streamoff>(originalLineLength + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        if (tokens.size() >= 6) {
            try {
                std::string type1 = tokens[0];
                std::string type2 = tokens[1];
                std::string type3 = tokens[2];
                std::string type4 = tokens[3];
                double kpsi = PrmParserStructures::safe_stod(tokens[4], "improper Kpsi");
                double psi0 = PrmParserStructures::safe_stod(tokens[5], "improper psi0");
                
                ff.add_improper_params(type1, type2, type3, type4, kpsi, psi0);
            } catch (const std::exception& e) {
                if (PrmParserSections::getDebugFlag()) std::cerr << "Warning: Skipping improper line due to parsing error: " << line << std::endl;
                continue;
            }
        }
    }
}

void PrmParserSections::parseNBFixSection(std::istream& input, pygcmc::model::ForceField& ff) {
    std::string line;
    while (std::getline(input, line)) {
        if (PrmParserStructures::isCommentLine(line)) continue;
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) continue;
        
        // Skip topology lines from STR files
        if (PrmParserStructures::isTopologyLine(line)) {
            if (PrmParserSections::getDebugFlag()) std::cerr << "Skipping topology line in NBFIX section: " << line << std::endl;
            continue;
        }
        
        // Check for section end or new section
        if (line == "END" || line == "end" || PrmParserStructures::isAtomsSection(line) || PrmParserStructures::isBondsSection(line) || 
            PrmParserStructures::isAnglesSection(line) || PrmParserStructures::isDihedralsSection(line) || 
            PrmParserStructures::isImproperSection(line) || PrmParserStructures::isNonbondedSection(line) || 
            line == "BOMLEV" || line == "WRNLEV" || line == "return") {
            input.seekg(-static_cast<std::streamoff>(line.length() + 1), std::ios::cur);
            break;
        }
        
        auto tokens = PrmParserStructures::tokenize(line);
        // Skip special directive lines like "HBOND CUTHB 0.5"
        if (tokens.size() >= 1 && (tokens[0] == "HBOND" || tokens[0] == "NBFIX")) {
            continue;
        }
        
        if (tokens.size() >= 4) {
            try {
                std::string type1 = tokens[0];
                std::string type2 = tokens[1];
                double epsilon = PrmParserStructures::safe_stod(tokens[2], "NBFIX epsilon for " + type1 + "-" + type2);
                double rmin = PrmParserStructures::safe_stod(tokens[3], "NBFIX Rmin for " + type1 + "-" + type2);
                
                // In CHARMM, NBFIX parameters are specified with full Rmin value
                // No need to multiply by 2 since we store the full Rmin value
                ff.add_nbfix(type1, type2, epsilon, rmin);
                
                if (PrmParserSections::getDebugFlag()) std::cerr << "Stored NBFIX for " << type1 << "-" << type2 
                    << ": epsilon = " << epsilon << ", Rmin = " << rmin << std::endl;
            } catch (const std::exception& e) {
                if (PrmParserSections::getDebugFlag()) std::cerr << "Warning: Skipping NBFIX line due to parsing error: " << line << std::endl;
            }
        }
    }
}

void PrmParserSections::parseNonbondedSection(std::istream& input, pygcmc::model::ForceField& ff, const std::string& firstLine) {
    std::string line = firstLine;
    std::string fullLine = PrmParserStructures::readContinuationLine(input, line);
    
    if (PrmParserSections::getDebugFlag()) std::cerr << "=== Entering NONBONDED section parsing ===" << std::endl;
    if (PrmParserSections::getDebugFlag()) std::cerr << "First line: [" << firstLine << "]" << std::endl;
    
    // Parse header parameters
    auto tokens = PrmParserStructures::tokenize(fullLine);
    if (!tokens.empty()) {
        if (PrmParserSections::getDebugFlag()) {
            std::cerr << "Initial tokens:";
            for (const auto& token : tokens) {
                std::cerr << " [" << token << "]";
            }
            std::cerr << std::endl;
        }
        
        // Skip the NONBONDED keyword
        if (tokens[0] == "NONBONDED") {
            if (PrmParserSections::getDebugFlag()) std::cerr << "Skipping NONBONDED keyword" << std::endl;
            tokens.erase(tokens.begin());
        }
        
        // Process parameters
        auto& params = ff.get_nonbonded_params();
        for (size_t i = 0; i < tokens.size(); ++i) {
            if (PrmParserSections::getDebugFlag()) std::cerr << "Processing parameter token: [" << tokens[i] << "]" << std::endl;
            
            if (tokens[i] == "nbxmod" && i + 1 < tokens.size()) {
                params.nbxmod = PrmParserStructures::safe_stoi(tokens[++i], "nbxmod");
                if (PrmParserSections::getDebugFlag()) std::cerr << "Set nbxmod = " << params.nbxmod << std::endl;
            } else if (tokens[i] == "cutnb" && i + 1 < tokens.size()) {
                params.cutnb = PrmParserStructures::safe_stod(tokens[++i], "cutnb");
                if (PrmParserSections::getDebugFlag()) std::cerr << "Set cutnb = " << params.cutnb << std::endl;
            } else if (tokens[i] == "ctofnb" && i + 1 < tokens.size()) {
                params.ctofnb = PrmParserStructures::safe_stod(tokens[++i], "ctofnb");
                if (PrmParserSections::getDebugFlag()) std::cerr << "Set ctofnb = " << params.ctofnb << std::endl;
            } else if (tokens[i] == "ctonnb" && i + 1 < tokens.size()) {
                params.ctonnb = PrmParserStructures::safe_stod(tokens[++i], "ctonnb");
                if (PrmParserSections::getDebugFlag()) std::cerr << "Set ctonnb = " << params.ctonnb << std::endl;
            } else if (tokens[i] == "eps" && i + 1 < tokens.size()) {
                params.eps = PrmParserStructures::safe_stod(tokens[++i], "eps");
                if (PrmParserSections::getDebugFlag()) std::cerr << "Set eps = " << params.eps << std::endl;
            } else if (tokens[i] == "e14fac" && i + 1 < tokens.size()) {
                params.e14fac = PrmParserStructures::safe_stod(tokens[++i], "e14fac");
                if (PrmParserSections::getDebugFlag()) std::cerr << "Set e14fac = " << params.e14fac << std::endl;
            } else if (tokens[i] == "wmin" && i + 1 < tokens.size()) {
                params.wmin = PrmParserStructures::safe_stod(tokens[++i], "wmin");
                if (PrmParserSections::getDebugFlag()) std::cerr << "Set wmin = " << params.wmin << std::endl;
            } else if (tokens[i] == "cdiel") {
                params.cdiel = true;
                if (PrmParserSections::getDebugFlag()) std::cerr << "Set cdiel = true" << std::endl;
            } else if (tokens[i] == "fshift") {
                params.fshift = true;
                if (PrmParserSections::getDebugFlag()) std::cerr << "Set fshift = true" << std::endl;
            } else if (tokens[i] == "vatom") {
                params.vatom = true;
                if (PrmParserSections::getDebugFlag()) std::cerr << "Set vatom = true" << std::endl;
            } else if (tokens[i] == "vdistance") {
                params.vdistance = true;
                if (PrmParserSections::getDebugFlag()) std::cerr << "Set vdistance = true" << std::endl;
            } else if (tokens[i] == "vfswitch") {
                params.vfswitch = true;
                if (PrmParserSections::getDebugFlag()) std::cerr << "Set vfswitch = true" << std::endl;
            }
        }
    }
    
    if (PrmParserSections::getDebugFlag()) std::cerr << "\n=== Starting atom type parameters parsing ===" << std::endl;
    if (PrmParserSections::getDebugFlag()) std::cerr << "Current lj_params map size: " << ff.get_num_lj_params() << std::endl;
    
    // Parse atom type parameters
    while (std::getline(input, line)) {
        if (PrmParserStructures::isCommentLine(line)) continue;
        
        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) {
            if (PrmParserSections::getDebugFlag()) std::cerr << "Skipping empty line" << std::endl;
            continue;
        }
        
        // Skip topology lines from STR files
        if (PrmParserStructures::isTopologyLine(line)) {
            if (PrmParserSections::getDebugFlag()) std::cerr << "Skipping topology line in nonbonded section: " << line << std::endl;
            continue;
        }
        
        if (PrmParserSections::getDebugFlag()) std::cerr << "Processing cleaned line: [" << line << "]" << std::endl;
        tokens = PrmParserStructures::tokenize(line);
        if (PrmParserSections::getDebugFlag()) {
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
            if (PrmParserSections::getDebugFlag()) std::cerr << "Found section end marker: " << tokens[0] << std::endl;
            break;
        }
        
        // If we find NBFIX, parse it as a new section
        if (PrmParserStructures::isNBFixSection(line)) {
            if (PrmParserSections::getDebugFlag()) std::cerr << "Found NBFIX section" << std::endl;
            parseNBFixSection(input, ff);
            break;
        }
        
        // Format: atomType ignored epsilon Rmin [ignored ignored ignored ignored]
        if (tokens.size() >= 4) {
            std::string atomType = tokens[0];
            // Skip the ignored value (usually 0.0)
            try {
                double epsilon = PrmParserStructures::safe_stod(tokens[2], "LJ epsilon for " + atomType);
                double rmin_half = PrmParserStructures::safe_stod(tokens[3], "LJ Rmin/2 for " + atomType);
                
                if (PrmParserSections::getDebugFlag()) {
                    std::cerr << "\n*** Parsing atom type: " << atomType << " ***" << std::endl;
                    std::cerr << "  epsilon = " << epsilon << std::endl;
                    std::cerr << "  rmin_half = " << rmin_half << std::endl;
                }
                
                ff.add_lj_params(atomType, epsilon, rmin_half);
                
                if (PrmParserSections::getDebugFlag()) {
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
    
    if (PrmParserSections::getDebugFlag()) std::cerr << "\n=== Finished NONBONDED section parsing ===" << std::endl;
    if (PrmParserSections::getDebugFlag()) std::cerr << "Final lj_params map size: " << ff.get_num_lj_params() << std::endl;
}

} // namespace io
} // namespace pygcmc