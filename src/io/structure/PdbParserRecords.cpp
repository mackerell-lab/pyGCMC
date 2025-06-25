// src/io/structure/PdbParserRecords.cpp

#include "PdbParserRecords.hpp"
#include <sstream>
#include <stdexcept>
#include <iostream>

namespace pygcmc {
namespace io {
namespace structure {

bool PdbParserRecords::parseAtomRecord(const std::string& line, 
                                      PdbParserStructures::RecordType type,
                                      model::Structure& structure,
                                      std::shared_ptr<model::Residue>& currentResidue) {
    try {
        if (line.length() < 54) return false;
        
        // Parse atom fields according to PDB format
        int serialNum = std::stoi(line.substr(6, 5));
        
        // Extract and trim atom name
        std::string atomName = line.substr(12, 4);
        size_t start = atomName.find_first_not_of(" ");
        size_t end = atomName.find_last_not_of(" ");
        if (start != std::string::npos && end != std::string::npos) {
            atomName = atomName.substr(start, end - start + 1);
        } else {
            atomName = "";
        }

        // Extract and trim residue name
        std::string resName = line.substr(17, 4);
        start = resName.find_first_not_of(" ");
        end = resName.find_last_not_of(" ");
        if (start != std::string::npos && end != std::string::npos) {
            resName = resName.substr(start, end - start + 1);
        } else {
            resName = "";
        }

        // Extract other fields
        std::string altLoc = line.length() > 16 ? std::string(1, line[16]) : "";
        
        std::string chainId = line.length() > 21 ? std::string(1, line[21]) : "";
        // Don't strip chain IDs - they can be single characters including space
        // The original logic handles empty chains correctly in char conversion below
        
        int resSeq = std::stoi(line.substr(22, 4));
        
        std::string iCode = line.length() > 26 ? std::string(1, line[26]) : "";
        if (!iCode.empty()) {
            size_t start = iCode.find_first_not_of(" ");
            if (start != std::string::npos) {
                iCode = iCode.substr(start);
            } else {
                iCode = "";
            }
        }
        
        double x = std::stod(line.substr(30, 8));
        double y = std::stod(line.substr(38, 8));
        double z = std::stod(line.substr(46, 8));
        
        double occupancy = (line.length() > 59) ? std::stod(line.substr(54, 6)) : 1.0;
        double tempFactor = (line.length() > 65) ? std::stod(line.substr(60, 6)) : 0.0;
        
        std::string element = (line.length() > 77) ? line.substr(76, 2) : "";
        start = element.find_first_not_of(" ");
        end = element.find_last_not_of(" ");
        if (start != std::string::npos && end != std::string::npos) {
            element = element.substr(start, end - start + 1);
        } else {
            element = "";
        }

        // Create residue if needed
        char chainChar = chainId.empty() ? ' ' : chainId[0];
        char iCodeChar = iCode.empty() ? ' ' : iCode[0];
        if (!currentResidue || currentResidue->get_resname() != resName || 
            currentResidue->get_chain() != chainChar || 
            currentResidue->get_ires() != resSeq || 
            currentResidue->get_inscode() != iCodeChar) {
            
            if (currentResidue) {
                currentResidue->calculate_center_of_mass();
            }
            
            currentResidue = std::make_shared<model::Residue>(resName, resSeq, "", 0, chainChar, iCodeChar);
            currentResidue->set_hetatm(type == PdbParserStructures::RecordType::HETATM);
            structure.add_residue(currentResidue);
        }

        // Create atom with default constructor and set properties (like original)
        auto atom = std::make_shared<model::Atom>();
        atom->set_bynu(serialNum);
        atom->set_type(atomName);
        atom->set_altloc(altLoc.empty() ? ' ' : altLoc[0]);
        atom->set_resname(resName);
        atom->set_chain(chainChar);
        atom->set_ires(resSeq);
        atom->set_inscode(iCodeChar);
        atom->set_coor(x, y, z);
        atom->set_occupancy(occupancy);
        atom->set_tempfactor(tempFactor);
        
        // Set HETATM flag if needed
        if (type == PdbParserStructures::RecordType::HETATM) {
            atom->set_hetatm(true);
        }
        
        // Set mass and charge from element masses
        const auto& masses = PdbParserStructures::getElementMasses();
        double mass = 0.0;
        auto it = masses.find(element);
        if (it != masses.end()) {
            mass = it->second;
        } else if (!atomName.empty()) {
            // Guess element from atom name
            std::string guessedElement = atomName.substr(0, 1);
            auto it2 = masses.find(guessedElement);
            if (it2 != masses.end()) {
                mass = it2->second;
            }
        }
        atom->set_mass_charge(mass, 0.0);

        // Add atom to BOTH structure and residue (like original parser)
        structure.add_atom(atom);
        currentResidue->add_atom(atom);
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error parsing ATOM record: " << e.what() << std::endl;
        return false;
    }
}

bool PdbParserRecords::parseTerRecord(const std::string& line,
                                     std::shared_ptr<model::Residue>& currentResidue,
                                     model::Structure& structure) {
    (void)line; (void)structure; // Suppress unused parameter warnings
    try {
        if (currentResidue) {
            currentResidue->calculate_center_of_mass();
            currentResidue = nullptr;
        }
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing TER record: " << e.what() << std::endl;
        return false;
    }
}

bool PdbParserRecords::parseHelixRecord(const std::string& line, model::Structure& structure) {
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

bool PdbParserRecords::parseSheetRecord(const std::string& line, model::Structure& structure) {
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

bool PdbParserRecords::parseSSBondRecord(const std::string& line, model::Structure& structure) {
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

bool PdbParserRecords::parseCryst1Record(const std::string& line, model::Structure& structure) {
    try {
        if (line.length() < 54) return false;
        
        double a = std::stod(line.substr(6, 9));
        double b = std::stod(line.substr(15, 9));
        double c = std::stod(line.substr(24, 9));
        double alpha = std::stod(line.substr(33, 7));
        double beta = std::stod(line.substr(40, 7));
        double gamma = std::stod(line.substr(47, 7));
        
        // Set box dimensions - include all 6 values
        structure.set_box_dimensions({a, b, c, alpha, beta, gamma});
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error parsing CRYST1 record: " << e.what() << std::endl;
        return false;
    }
}

} // namespace structure
} // namespace io
} // namespace pygcmc