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
                parseTerRecord(line, currentResidue, result);
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
        
        // Extract atom name and trim leading/trailing spaces
        std::string atomName = line.substr(12, 4);
        if (!atomName.empty()) {
            size_t start = atomName.find_first_not_of(" ");
            size_t end = atomName.find_last_not_of(" ");
            if (start != std::string::npos && end != std::string::npos) {
                atomName = atomName.substr(start, end - start + 1);
            } else {
                atomName = "";
            }
        }

        // Extract residue name and trim leading/trailing spaces
        std::string resName = line.substr(17, 4);
        if (!resName.empty()) {
            size_t start = resName.find_first_not_of(" ");
            size_t end = resName.find_last_not_of(" ");
            if (start != std::string::npos && end != std::string::npos) {
                resName = resName.substr(start, end - start + 1);
            } else {
                resName = "";
            }
        }

        // Extract and trim other fields
        std::string altLoc = line.length() > 16 ? std::string(1, line[16]) : "";
        if (!altLoc.empty()) {
            size_t start = altLoc.find_first_not_of(" ");
            if (start != std::string::npos) {
                altLoc = altLoc.substr(start);
            }
        }
        
        std::string chainId = line.length() > 21 ? std::string(1, line[21]) : "";
        if (!chainId.empty()) {
            size_t start = chainId.find_first_not_of(" ");
            if (start != std::string::npos) {
                chainId = chainId.substr(start);
            }
        }
        
        int resSeq = std::stoi(line.substr(22, 4));
        
        std::string iCode = line.length() > 26 ? std::string(1, line[26]) : "";
        if (!iCode.empty()) {
            size_t start = iCode.find_first_not_of(" ");
            if (start != std::string::npos) {
                iCode = iCode.substr(start);
            }
        }
        
        double x = std::stod(line.substr(30, 8));
        double y = std::stod(line.substr(38, 8));
        double z = std::stod(line.substr(46, 8));
        
        double occupancy = (line.length() > 59) ? 
            std::stod(line.substr(54, 6)) : 1.0;
        double tempFactor = (line.length() > 65) ? 
            std::stod(line.substr(60, 6)) : 0.0;
        
        std::string segId = (line.length() > 75) ? 
            line.substr(72, 4) : "";
        if (!segId.empty()) {
            size_t start = segId.find_first_not_of(" ");
            size_t end = segId.find_last_not_of(" ");
            if (start != std::string::npos && end != std::string::npos) {
                segId = segId.substr(start, end - start + 1);
            } else {
                segId = "";
            }
        }
        
        std::string element = (line.length() > 77) ? 
            line.substr(76, 2) : "";
        if (!element.empty()) {
            size_t start = element.find_first_not_of(" ");
            size_t end = element.find_last_not_of(" ");
            if (start != std::string::npos && end != std::string::npos) {
                element = element.substr(start, end - start + 1);
            } else {
                element = "";
            }
        }
        
        std::string charge = (line.length() > 79) ? 
            line.substr(78, 2) : "";
        if (!charge.empty()) {
            size_t start = charge.find_first_not_of(" ");
            size_t end = charge.find_last_not_of(" ");
            if (start != std::string::npos && end != std::string::npos) {
                charge = charge.substr(start, end - start + 1);
            } else {
                charge = "";
            }
        }

        // Create atom with stripped values
        auto atom = std::make_shared<model::Atom>();
        atom->setBynu(serialNum);
        atom->setType(atomName);
        atom->setAltloc(altLoc.empty() ? ' ' : altLoc[0]);
        atom->setResname(resName);
        atom->setChain(chainId.empty() ? ' ' : chainId[0]);
        atom->setIres(resSeq);
        atom->setInscode(iCode.empty() ? ' ' : iCode[0]);
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
        if (!currentResidue || 
            currentResidue->getChain() != chainId[0] ||
            currentResidue->getIres() != resSeq ||
            currentResidue->getInscode() != iCode[0]) {
            // Check if residue already exists
            bool found = false;
            for (auto& res : result.residues) {
                if (res->getChain() == chainId[0] &&
                    res->getIres() == resSeq &&
                    res->getInscode() == iCode[0]) {
                    currentResidue = res;
                    found = true;
                    break;
                }
            }
            if (!found) {
                // Create new residue with current name and HETATM status
                currentResidue = std::make_shared<model::Residue>(
                    resName, resSeq, segId, 0, chainId[0], iCode[0]);
                currentResidue->setHetatm(type == RecordType::HETATM);
                result.residues.push_back(currentResidue);
            } else if (currentResidue->getResname() != resName || 
                      currentResidue->isHetatm() != (type == RecordType::HETATM)) {
                // If residue exists but has different name/HETATM status,
                // create a new one with updated properties
                auto newResidue = std::make_shared<model::Residue>(
                    resName, resSeq, segId, 0, chainId[0], iCode[0]);
                newResidue->setHetatm(type == RecordType::HETATM);
                // Copy existing atoms
                for (const auto& existingAtom : currentResidue->getAtoms()) {
                    newResidue->addAtom(existingAtom);
                }
                // Replace old residue with new one
                auto it = std::find(result.residues.begin(), result.residues.end(), currentResidue);
                if (it != result.residues.end()) {
                    *it = newResidue;
                }
                currentResidue = newResidue;
            }
        }
        currentResidue->addAtom(atom);

    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing ATOM/HETATM record: " + 
                               std::string(e.what()));
    }
}

void PDBParser::parseTerRecord(const std::string& line,
                             std::shared_ptr<model::Residue>& currentResidue,
                             ParseResult& result) {
    try {
        std::string chainId = line.length() > 21 ? std::string(1, line[21]) : "";
        if (!chainId.empty()) {
            size_t start = chainId.find_first_not_of(" ");
            if (start != std::string::npos) {
                chainId = chainId.substr(start);
            }
        }
        
        int resSeq = 0;
        std::string iCode = " ";
        std::string resName;

        // Get residue information from current residue if available
        if (currentResidue) {
            chainId = std::string(1, currentResidue->getChain());
            resSeq = currentResidue->getIres();
            iCode = std::string(1, currentResidue->getInscode());
            resName = currentResidue->getResname();
        }

        // Try to parse residue sequence number if present
        if (line.length() >= 26) {
            try {
                resSeq = std::stoi(line.substr(22, 4));
            } catch (const std::exception&) {
                // Use default or current residue value
            }
        }

        // Store terminal information
        result.terminals.push_back(TerminalInfo{
            chainId.empty() ? ' ' : chainId[0], 
            resSeq, 
            iCode.empty() ? ' ' : iCode[0], 
            resName
        });

        // Reset current residue pointer
        currentResidue = nullptr;
    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing TER record: " + std::string(e.what()));
    }
}

void PDBParser::parseHelixRecord(const std::string& line, ParseResult& result) {
    try {
        std::string helixId = line.substr(11, 3);
        if (!helixId.empty()) {
            size_t start = helixId.find_first_not_of(" ");
            size_t end = helixId.find_last_not_of(" ");
            if (start != std::string::npos && end != std::string::npos) {
                helixId = helixId.substr(start, end - start + 1);
            } else {
                helixId = "";
            }
        }
        
        std::string initResName = line.substr(15, 3);
        if (!initResName.empty()) {
            size_t start = initResName.find_first_not_of(" ");
            size_t end = initResName.find_last_not_of(" ");
            if (start != std::string::npos && end != std::string::npos) {
                initResName = initResName.substr(start, end - start + 1);
            } else {
                initResName = "";
            }
        }
        
        std::string initChainId = std::string(1, line[19]);
        if (!initChainId.empty()) {
            size_t start = initChainId.find_first_not_of(" ");
            if (start != std::string::npos) {
                initChainId = initChainId.substr(start);
            }
        }
        
        std::string initICode = line.length() > 25 ? std::string(1, line[25]) : "";
        if (!initICode.empty()) {
            size_t start = initICode.find_first_not_of(" ");
            if (start != std::string::npos) {
                initICode = initICode.substr(start);
            }
        }
        
        std::string endResName = line.substr(27, 3);
        if (!endResName.empty()) {
            size_t start = endResName.find_first_not_of(" ");
            size_t end = endResName.find_last_not_of(" ");
            if (start != std::string::npos && end != std::string::npos) {
                endResName = endResName.substr(start, end - start + 1);
            } else {
                endResName = "";
            }
        }
        
        std::string endChainId = std::string(1, line[31]);
        if (!endChainId.empty()) {
            size_t start = endChainId.find_first_not_of(" ");
            if (start != std::string::npos) {
                endChainId = endChainId.substr(start);
            }
        }
        
        std::string endICode = line.length() > 37 ? std::string(1, line[37]) : "";
        if (!endICode.empty()) {
            size_t start = endICode.find_first_not_of(" ");
            if (start != std::string::npos) {
                endICode = endICode.substr(start);
            }
        }
        
        int helixClass = std::stoi(line.substr(38, 2));
        
        // Store helix information
        result.helices[initChainId].push_back(helixClass);

    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing HELIX record: " + std::string(e.what()));
    }
}

void PDBParser::parseSheetRecord(const std::string& line, ParseResult& result) {
    try {
        if (line.length() < 38) {
            throw std::runtime_error("SHEET record too short");
        }

        std::string strandStr = line.substr(7, 3);
        if (!strandStr.empty()) {
            size_t start = strandStr.find_first_not_of(" ");
            size_t end = strandStr.find_last_not_of(" ");
            if (start != std::string::npos && end != std::string::npos) {
                strandStr = strandStr.substr(start, end - start + 1);
            } else {
                strandStr = "";
            }
        }
        
        std::string sheetId = line.substr(11, 3);
        if (!sheetId.empty()) {
            size_t start = sheetId.find_first_not_of(" ");
            size_t end = sheetId.find_last_not_of(" ");
            if (start != std::string::npos && end != std::string::npos) {
                sheetId = sheetId.substr(start, end - start + 1);
            } else {
                sheetId = "";
            }
        }
        
        std::string initChainId = std::string(1, line[21]);
        if (!initChainId.empty()) {
            size_t start = initChainId.find_first_not_of(" ");
            if (start != std::string::npos) {
                initChainId = initChainId.substr(start);
            }
        }
        
        int strandNum = strandStr.empty() ? 1 : std::stoi(strandStr);
        
        int sense = 0;
        if (line.length() >= 40) {
            std::string senseStr = line.substr(38, 2);
            if (!senseStr.empty()) {
                size_t start = senseStr.find_first_not_of(" ");
                size_t end = senseStr.find_last_not_of(" ");
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
        }
        
        std::string sheetInfo = sheetId + ":" + 
            std::to_string(strandNum) + ":" + 
            std::to_string(sense);
        
        if (result.sheets.find(initChainId) == result.sheets.end()) {
            result.sheets[initChainId] = std::vector<std::string>();
        }
        result.sheets[initChainId].push_back(sheetInfo);

    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing SHEET record: " + std::string(e.what()));
    }
}

void PDBParser::parseSSBondRecord(const std::string& line, ParseResult& result) {
    try {
        std::string chain1 = std::string(1, line[15]);
        if (!chain1.empty()) {
            size_t start = chain1.find_first_not_of(" ");
            if (start != std::string::npos) {
                chain1 = chain1.substr(start);
            }
        }
        
        int resnum1 = std::stoi(line.substr(17, 4));
        
        std::string inscode1 = line.length() > 21 ? std::string(1, line[21]) : "";
        if (!inscode1.empty()) {
            size_t start = inscode1.find_first_not_of(" ");
            if (start != std::string::npos) {
                inscode1 = inscode1.substr(start);
            }
        }
        
        std::string chain2 = std::string(1, line[29]);
        if (!chain2.empty()) {
            size_t start = chain2.find_first_not_of(" ");
            if (start != std::string::npos) {
                chain2 = chain2.substr(start);
            }
        }
        
        int resnum2 = std::stoi(line.substr(31, 4));
        
        std::string inscode2 = line.length() > 35 ? std::string(1, line[35]) : "";
        if (!inscode2.empty()) {
            size_t start = inscode2.find_first_not_of(" ");
            if (start != std::string::npos) {
                inscode2 = inscode2.substr(start);
            }
        }
        
        std::string bondInfo = chain1 + ":" + 
            std::to_string(resnum1) + inscode1 + "-" +
            chain2 + ":" + 
            std::to_string(resnum2) + inscode2;
        
        result.ssbonds.push_back(bondInfo);

    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing SSBOND record: " + std::string(e.what()));
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