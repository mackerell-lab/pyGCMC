// src/io/pdbParser.cpp

#include "pdbParser.hpp"
#include "model/atom.hpp"
#include "model/residue.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace pygcmc {
namespace io {

PDBParser::ParseResult PDBParser::parseFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Could not open PDB file: " + filename);
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return parseString(buffer.str());
}

PDBParser::ParseResult PDBParser::parseString(const std::string& pdbStr) {
    ParseResult result;
    std::istringstream iss(pdbStr);
    std::string line;
    
    // Current chain tracking
    std::string currentChain;
    std::shared_ptr<model::Residue> currentResidue;
    
    while (std::getline(iss, line)) {
        if (line.length() < 6) continue;
        
        RecordType recordType = getRecordType(line);
        switch (recordType) {
            case RecordType::ATOM:
            case RecordType::HETATM:
                parseAtomRecord(line, recordType, result, currentResidue);
                break;
            case RecordType::TER:
                parseTerRecord(line, currentResidue);
                break;
            case RecordType::HELIX:
                parseHelixRecord(line, result);
                break;
            case RecordType::SHEET:
                parseSheetRecord(line, result);
                break;
            case RecordType::SSBOND:
                parseSSBondRecord(line, result);
                break;
            case RecordType::CRYST1:
                parseCryst1Record(line, result);
                break;
            default:
                continue;
        }
    }
    
    return result;
}

PDBParser::RecordType PDBParser::getRecordType(const std::string& line) {
    std::string recordName = line.substr(0, 6);
    if (recordName == "ATOM  ") return RecordType::ATOM;
    if (recordName == "HETATM") return RecordType::HETATM;
    if (recordName.substr(0, 3) == "TER") return RecordType::TER;
    if (recordName.substr(0, 5) == "HELIX") return RecordType::HELIX;
    if (recordName.substr(0, 5) == "SHEET") return RecordType::SHEET;
    if (recordName.substr(0, 6) == "SSBOND") return RecordType::SSBOND;
    if (recordName.substr(0, 6) == "CRYST1") return RecordType::CRYST1;
    return RecordType::UNKNOWN;
}

void PDBParser::parseAtomRecord(const std::string& line, RecordType type,
                              ParseResult& result,
                              std::shared_ptr<model::Residue>& currentResidue) {
    try {
        // Parse atom fields according to PDB format
        int serialNum = std::stoi(line.substr(6, 5));
        std::string atomName = line.substr(12, 4);
        char altLoc = (line.length() > 16) ? line[16] : ' ';
        std::string resName = line.substr(17, 4);
        // Trim trailing whitespace from residue name
        resName.erase(std::find_if(resName.rbegin(), resName.rend(), 
            [](unsigned char ch) { return !std::isspace(ch); }).base(), resName.end());
        char chainId = (line.length() > 21) ? line[21] : ' ';
        int resSeq = std::stoi(line.substr(22, 4));
        char iCode = (line.length() > 26) ? line[26] : ' ';
        
        double x = std::stod(line.substr(30, 8));
        double y = std::stod(line.substr(38, 8));
        double z = std::stod(line.substr(46, 8));
        
        double occupancy = (line.length() > 59) ? 
            std::stod(line.substr(54, 6)) : 1.0;
        double tempFactor = (line.length() > 65) ? 
            std::stod(line.substr(60, 6)) : 0.0;
        
        std::string segId = (line.length() > 75) ? 
            line.substr(72, 4) : "";
        std::string element = (line.length() > 77) ? 
            line.substr(76, 2) : "";
        std::string charge = (line.length() > 79) ? 
            line.substr(78, 2) : "";

        // Only use record type to determine HETATM status
        bool isHetatm = type == RecordType::HETATM;

        // Create atom
        auto atom = std::make_shared<model::Atom>();
        atom->setBynu(serialNum);
        atom->setType(atomName);
        atom->setAltloc(altLoc);
        atom->setResname(resName);
        atom->setChain(chainId);
        atom->setIres(resSeq);
        atom->setInscode(iCode);
        atom->setCoor(x, y, z);
        atom->setOccupancy(occupancy);
        atom->setTempfactor(tempFactor);
        atom->setSegid(segId);
        atom->setElement(element);
        atom->setChargeString(charge);
        atom->setHetatm(isHetatm);

        // Add to result
        result.atoms.push_back(atom);

        // Handle residue
        if (!currentResidue || 
            currentResidue->getChain() != chainId ||
            currentResidue->getIres() != resSeq ||
            currentResidue->getInscode() != iCode ||
            currentResidue->getResname() != resName ||
            currentResidue->isHetatm() != isHetatm) {
            // Check if residue already exists
            bool found = false;
            for (auto& res : result.residues) {
                if (res->getChain() == chainId &&
                    res->getIres() == resSeq &&
                    res->getInscode() == iCode &&
                    res->getResname() == resName &&
                    res->isHetatm() == isHetatm) {
                    currentResidue = res;
                    found = true;
                    break;
                }
            }
            if (!found) {
                // Create new residue
                currentResidue = std::make_shared<model::Residue>(
                    resName, resSeq, segId, 0, chainId, iCode);
                currentResidue->setHetatm(isHetatm);
                result.residues.push_back(currentResidue);
            }
        }
        currentResidue->addAtom(atom);

    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing ATOM/HETATM record: " + 
                               std::string(e.what()));
    }
}

void PDBParser::parseTerRecord([[maybe_unused]] const std::string& line,
                              std::shared_ptr<model::Residue>& currentResidue) {
    // Reset current residue pointer to indicate end of chain
    currentResidue = nullptr;
}

void PDBParser::parseHelixRecord(const std::string& line, ParseResult& result) {
    try {
        [[maybe_unused]] int serialNum = std::stoi(line.substr(7, 3));
        std::string helixId = line.substr(11, 3);
        std::string initResName = line.substr(15, 3);
        char initChainId = line[19];
        [[maybe_unused]] int initSeqNum = std::stoi(line.substr(21, 4));
        [[maybe_unused]] char initICode = line[25];
        std::string endResName = line.substr(27, 3);
        [[maybe_unused]] char endChainId = line[31];
        [[maybe_unused]] int endSeqNum = std::stoi(line.substr(33, 4));
        [[maybe_unused]] char endICode = line[37];
        int helixClass = std::stoi(line.substr(38, 2));
        
        // Store helix information
        result.helices[std::string(1, initChainId)].push_back(helixClass);

    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing HELIX record: " + 
                               std::string(e.what()));
    }
}

void PDBParser::parseSheetRecord(const std::string& line, ParseResult& result) {
    try {
        // Ensure line is long enough for minimal SHEET record
        if (line.length() < 38) {
            throw std::runtime_error("SHEET record too short");
        }

        // Parse fields according to PDB format specification
        std::string strandStr = line.substr(7, 3);      // Strand number (8-10)
        std::string sheetId = line.substr(11, 3);       // Sheet ID (12-14)
        char initChainId = line[21];                    // Initial chain ID (22)
        
        // Trim whitespace
        strandStr.erase(std::remove_if(strandStr.begin(), strandStr.end(), ::isspace), strandStr.end());
        sheetId.erase(std::remove_if(sheetId.begin(), sheetId.end(), ::isspace), sheetId.end());
        
        // Convert to integer with validation
        int strandNum = strandStr.empty() ? 1 : std::stoi(strandStr);
        
        // Parse sense value (columns 39-40)
        int sense = 0;  // default for first strand
        if (line.length() >= 40) {
            std::string senseStr = line.substr(38, 2);  // columns 39-40
            // Print raw sense string for debugging
            std::cerr << "Raw sense string: '" << senseStr << "'" << std::endl;
            
            // Remove spaces but keep negative sign
            senseStr.erase(std::remove_if(senseStr.begin(), senseStr.end(), 
                [](unsigned char c) { return std::isspace(c); }), senseStr.end());
            
            // Print trimmed sense string
            std::cerr << "Trimmed sense string: '" << senseStr << "'" << std::endl;
            
            if (!senseStr.empty()) {
                try {
                    sense = std::stoi(senseStr);
                    std::cerr << "Converted sense value: " << sense << std::endl;
                } catch (const std::exception& e) {
                    std::cerr << "Failed to convert sense string: " << e.what() << std::endl;
                    // Keep default sense value if conversion fails
                }
            }
        }
        
        // Create sheet info string (format: sheetId:strandNum:sense)
        std::string sheetInfo = sheetId + ":" + 
            std::to_string(strandNum) + ":" + 
            std::to_string(sense);
        
        // Initialize vector if chain not present
        if (result.sheets.find(std::string(1, initChainId)) == result.sheets.end()) {
            result.sheets[std::string(1, initChainId)] = std::vector<std::string>();
        }
        result.sheets[std::string(1, initChainId)].push_back(sheetInfo);

    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing SHEET record: " + std::string(e.what()));
    }
}

void PDBParser::parseSSBondRecord(const std::string& line, ParseResult& result) {
    try {
        [[maybe_unused]] int serialNum = std::stoi(line.substr(7, 3));
        
        // First CYS
        char chain1 = line[15];
        int resnum1 = std::stoi(line.substr(17, 4));
        char inscode1 = line[21];
        
        // Second CYS
        char chain2 = line[29];
        int resnum2 = std::stoi(line.substr(31, 4));
        char inscode2 = line[35];
        
        // Create bond info string
        std::string bondInfo = std::string(1, chain1) + ":" + 
            std::to_string(resnum1) + std::string(1, inscode1) + "-" +
            std::string(1, chain2) + ":" + 
            std::to_string(resnum2) + std::string(1, inscode2);
        
        result.ssbonds.push_back(bondInfo);

    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing SSBOND record: " + 
                               std::string(e.what()));
    }
}

void PDBParser::parseCryst1Record(const std::string& line, ParseResult& result) {
    try {
        // Parse unit cell parameters according to PDB format
        double a = std::stod(line.substr(6, 9));
        double b = std::stod(line.substr(15, 9));
        double c = std::stod(line.substr(24, 9));
        double alpha = std::stod(line.substr(33, 7));
        double beta = std::stod(line.substr(40, 7));
        double gamma = std::stod(line.substr(47, 7));
        
        // Store box dimensions
        result.boxDimensions = std::vector<double>{a, b, c, alpha, beta, gamma};
        
    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing CRYST1 record: " + 
                               std::string(e.what()));
    }
}

} // namespace io
} // namespace pygcmc 