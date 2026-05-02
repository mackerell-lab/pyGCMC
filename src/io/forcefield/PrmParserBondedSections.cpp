// src/io/forcefield/PrmParserBondedSections.cpp

#include "PrmParserBondedSections.hpp"
#include "PrmParserStructures.hpp"
#include "model/ModelModule.hpp"

namespace pygcmc {
namespace io {

void PrmParserBondedSections::parseDihedralsSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output) {
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

void PrmParserBondedSections::parseImproperSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output) {
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
                if (debug_output) std::cerr << "Warning: Skipping improper line due to parsing error: " << line << std::endl;
                continue;
            }
        }
    }
}

void PrmParserBondedSections::parseNBFixSection(std::istream& input, pygcmc::model::ForceField& ff, bool& debug_output) {
    std::string line;
    while (std::getline(input, line)) {
        if (PrmParserStructures::isCommentLine(line)) continue;

        line = PrmParserStructures::removeComments(line);
        line = PrmParserStructures::trim(line);
        if (line.empty()) continue;

        // Skip topology lines from STR files
        if (PrmParserStructures::isTopologyLine(line)) {
            if (debug_output) std::cerr << "Skipping topology line in NBFIX section: " << line << std::endl;
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

                if (debug_output) std::cerr << "Stored NBFIX for " << type1 << "-" << type2
                    << ": epsilon = " << epsilon << ", Rmin = " << rmin << std::endl;
            } catch (const std::exception& e) {
                if (debug_output) std::cerr << "Warning: Skipping NBFIX line due to parsing error: " << line << std::endl;
            }
        }
    }
}

} // namespace io
} // namespace pygcmc
