// src/io/structure/PdbParserStructureRecords.cpp

#include "PdbParserStructureRecords.hpp"
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <set>

namespace pygcmc {
namespace io {
namespace structure {

bool PdbParserStructureRecords::parseHelixRecord(const std::string& line, model::Structure& structure) {
    try {
        if (line.length() < 40) return false;

        // Parse HELIX record according to PDB format
        // HELIX serial helixID initResName initChainID initSeqNum initICode endResName endChainID endSeqNum endICode helixClass comment length

        std::string helixId = line.substr(11, 3);
        // Trim helixId
        size_t start = helixId.find_first_not_of(" ");
        size_t end = helixId.find_last_not_of(" ");
        if (start != std::string::npos && end != std::string::npos) {
            helixId = helixId.substr(start, end - start + 1);
        }

        std::string initResName = line.substr(15, 3);
        // Trim initResName
        start = initResName.find_first_not_of(" ");
        end = initResName.find_last_not_of(" ");
        if (start != std::string::npos && end != std::string::npos) {
            initResName = initResName.substr(start, end - start + 1);
        }

        char initChainId = line.length() > 19 ? line[19] : ' ';
        int initSeqNum = std::stoi(line.substr(21, 4));
        char initICode = line.length() > 25 ? line[25] : ' ';

        std::string endResName = line.substr(27, 3);
        // Trim endResName
        start = endResName.find_first_not_of(" ");
        end = endResName.find_last_not_of(" ");
        if (start != std::string::npos && end != std::string::npos) {
            endResName = endResName.substr(start, end - start + 1);
        }

        char endChainId = line.length() > 31 ? line[31] : ' ';
        int endSeqNum = std::stoi(line.substr(33, 4));
        char endICode = line.length() > 37 ? line[37] : ' ';

        // Parse helix class (column 39-40)
        int helixClass = 1; // Default to alpha helix
        if (line.length() > 39) {
            try {
                helixClass = std::stoi(line.substr(38, 2));
            } catch (...) {
                helixClass = 1; // Default if parsing fails
            }
        }

        // Create SecondaryStructure object
        model::Structure::SecondaryStructure helix;
        helix.id = helixId;
        helix.initResName = initResName;
        helix.initChainId = initChainId;
        helix.initSeqNum = initSeqNum;
        helix.initICode = initICode;
        helix.endResName = endResName;
        helix.endChainId = endChainId;
        helix.endSeqNum = endSeqNum;
        helix.endICode = endICode;
        helix.structureClass = helixClass;

        // Add to structure using chain ID as key
        std::string chainIdStr(1, initChainId);
        structure.add_helix(chainIdStr, helix);

        return true;

    } catch (const std::exception& e) {
        std::cerr << "Error parsing HELIX record: " << e.what() << std::endl;
        return false;
    }
}

bool PdbParserStructureRecords::parseSheetRecord(const std::string& line, model::Structure& structure) {
    try {
        if (line.length() < 33) return false;

        // Parse SHEET record according to PDB format
        std::string strandStr = line.substr(7, 3);
        size_t start = strandStr.find_first_not_of(" ");
        size_t end = strandStr.find_last_not_of(" ");
        if (start != std::string::npos && end != std::string::npos) {
            strandStr = strandStr.substr(start, end - start + 1);
        }

        std::string sheetId = line.substr(11, 3);
        start = sheetId.find_first_not_of(" ");
        end = sheetId.find_last_not_of(" ");
        if (start != std::string::npos && end != std::string::npos) {
            sheetId = sheetId.substr(start, end - start + 1);
        }

        std::string chainId = line.length() > 21 ? std::string(1, line[21]) : "";

        int strandNum = strandStr.empty() ? 1 : std::stoi(strandStr);

        // Parse sense (column 39-40)
        int sense = 0;
        if (line.length() >= 40) {
            std::string senseStr = line.substr(38, 2);
            start = senseStr.find_first_not_of(" ");
            end = senseStr.find_last_not_of(" ");
            if (start != std::string::npos && end != std::string::npos) {
                senseStr = senseStr.substr(start, end - start + 1);
                if (!senseStr.empty()) {
                    try {
                        sense = std::stoi(senseStr);
                    } catch (const std::exception&) {
                        // Keep default sense value if conversion fails
                    }
                }
            }
        }

        // Format sheet info like original: "sheetId:strandNum:sense"
        std::string sheetInfo = sheetId + ":" + std::to_string(strandNum) + ":" + std::to_string(sense);
        structure.add_sheet(chainId, sheetInfo);

        return true;

    } catch (const std::exception& e) {
        std::cerr << "Error parsing SHEET record: " << e.what() << std::endl;
        return false;
    }
}

bool PdbParserStructureRecords::parseSSBondRecord(const std::string& line, model::Structure& structure) {
    try {
        if (line.length() < 30) return false;

        // Parse SSBOND record according to PDB format
        // SSBOND serNum resName1 chainID1 seqNum1 icode1 resName2 chainID2 seqNum2 icode2

        char cys1Chain = line.length() > 15 ? line[15] : ' ';
        int cys1ResSeq = std::stoi(line.substr(17, 4));
        char cys1ICode = line.length() > 21 ? line[21] : ' ';

        char cys2Chain = line.length() > 29 ? line[29] : ' ';
        int cys2ResSeq = std::stoi(line.substr(31, 4));
        char cys2ICode = line.length() > 35 ? line[35] : ' ';

        // Format: "chainId:resSeq icode-chainId:resSeq icode" (matching original implementation)
        std::string icode1 = (cys1ICode != ' ') ? std::string(1, cys1ICode) : " ";
        std::string icode2 = (cys2ICode != ' ') ? std::string(1, cys2ICode) : " ";
        std::string bond = std::string(1, cys1Chain) + ":" + std::to_string(cys1ResSeq) + icode1 + "-" +
                          std::string(1, cys2Chain) + ":" + std::to_string(cys2ResSeq) + icode2;

        structure.add_ssbond(bond);
        return true;

    } catch (const std::exception& e) {
        std::cerr << "Error parsing SSBOND record: " << e.what() << std::endl;
        return false;
    }
}

bool PdbParserStructureRecords::parseCryst1Record(const std::string& line, model::Structure& structure) {
    try {
        // Many PDBs follow the fixed-column CRYST1 format, but some minimal/legacy writers emit
        // whitespace-separated values. Prefer robust token parsing to avoid silent mis-parses.
        std::istringstream iss(line);
        std::string record;
        iss >> record;
        if (record != "CRYST1") return false;

        double a = 0.0, b = 0.0, c = 0.0;
        if (!(iss >> a >> b >> c)) {
            return false;
        }

        // Angles are optional in some minimal files; default to orthorhombic.
        double alpha = 90.0, beta = 90.0, gamma = 90.0;
        if (!(iss >> alpha >> beta >> gamma)) {
            alpha = 90.0;
            beta = 90.0;
            gamma = 90.0;
        }

        structure.set_box_dimensions({a, b, c, alpha, beta, gamma});
        return true;

    } catch (const std::exception& e) {
        std::cerr << "Error parsing CRYST1 record: " << e.what() << std::endl;
        return false;
    }
}

bool PdbParserStructureRecords::isSmallMolecule(const std::string& resName) {
    // Common GCMC molecules and small molecules that should use continuity checking
    static const std::set<std::string> smallMolecules = {
        // GCMC molecules from 4wp7 system (actual names)
        "ACE", "BEN", "DME", "FOR", "IMI", "MAM", "MEO", "PRP",
        // Legacy names for compatibility
        "ACEY", "BENX", "DMEE", "FORM", "IMIA", "MAMY", "MEOH", "PRPX",
        // Common water and ions
        "SOL", "HOH", "WAT", "TIP", "TIP3", "SPC", "SPCE",
        "NA", "CL", "K", "CA", "MG", "ZN", "FE", "CU",
        "NA+", "CL-", "K+", "CA2+", "MG2+", "ZN2+", "FE2+", "FE3+", "CU2+",
        // Common small molecules and ligands
        "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "NAD", "FAD", "FMN",
        "HEM", "CHL", "BCL", "PHE", "TYR", "TRP", "HIS", "ARG", "LYS",
        // Organic solvents
        "DMSO", "ACE", "MOH", "EOH", "CHX", "DCM", "TCE", "TOL", "BEN",
        // Gases
        "CO2", "O2", "N2", "H2", "CO", "NH3", "CH4", "H2S", "SO2"
    };

    return smallMolecules.find(resName) != smallMolecules.end();
}

} // namespace structure
} // namespace io
} // namespace pygcmc
