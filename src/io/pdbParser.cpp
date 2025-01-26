// src/io/pdbParser.cpp

#include "pdbParser.hpp"
#include "model/atom.hpp"
#include "model/residue.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

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
        std::string resName = line.substr(17, 3);
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
        atom->setHetatm(type == RecordType::HETATM);

        // Add to result
        result.atoms.push_back(atom);

        // Handle residue
        std::string residueKey = chainId + std::to_string(resSeq) + iCode;
        if (!currentResidue || 
            currentResidue->getChain() != chainId ||
            currentResidue->getIres() != resSeq ||
            currentResidue->getInscode() != iCode) {
            // Create new residue
            currentResidue = std::make_shared<model::Residue>(
                resName, resSeq, segId, 0, chainId, iCode);
            result.residues.push_back(currentResidue);
        }
        currentResidue->addAtom(atom);

    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing ATOM/HETATM record: " + 
                               std::string(e.what()));
    }
}

void PDBParser::parseTerRecord(const std::string& line,
                             std::shared_ptr<model::Residue>& currentResidue) {
    currentResidue.reset();  // End current residue
}

void PDBParser::parseHelixRecord(const std::string& line, ParseResult& result) {
    try {
        int serialNum = std::stoi(line.substr(7, 3));
        std::string helixId = line.substr(11, 3);
        std::string initResName = line.substr(15, 3);
        char initChainId = line[19];
        int initSeqNum = std::stoi(line.substr(21, 4));
        char initICode = line[25];
        std::string endResName = line.substr(27, 3);
        char endChainId = line[31];
        int endSeqNum = std::stoi(line.substr(33, 4));
        char endICode = line[37];
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
        int strandNum = std::stoi(line.substr(7, 3));
        std::string sheetId = line.substr(11, 3);
        int numStrands = std::stoi(line.substr(14, 2));
        
        // Initial residue
        std::string initResName = line.substr(17, 3);
        char initChainId = line[21];
        int initSeqNum = std::stoi(line.substr(22, 4));
        char initICode = line[26];
        
        // Terminal residue
        std::string endResName = line.substr(28, 3);
        char endChainId = line[32];
        int endSeqNum = std::stoi(line.substr(33, 4));
        char endICode = line[37];
        
        // Sense
        int sense = std::stoi(line.substr(38, 2));
        
        // Create sheet info string
        std::string sheetInfo = sheetId + ":" + 
            std::to_string(strandNum) + ":" + 
            std::to_string(sense);
        
        result.sheets[std::string(1, initChainId)].push_back(sheetInfo);

    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing SHEET record: " + 
                               std::string(e.what()));
    }
}

void PDBParser::parseSSBondRecord(const std::string& line, ParseResult& result) {
    try {
        int serialNum = std::stoi(line.substr(7, 3));
        
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

} // namespace io
} // namespace pygcmc 