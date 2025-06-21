// src/io/pdbParser.cpp

#include "pdbParser.hpp"
#include "model/atom.hpp"
#include "model/residue.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <filesystem>

namespace pygcmc {
namespace io {

// Define the element mass table
const std::map<std::string, double> PDBParser::ELEMENT_MASSES = {
    {"H", 1.008},   // Hydrogen
    {"He", 4.003},  // Helium
    {"Li", 6.941},  // Lithium
    {"Be", 9.012},  // Beryllium
    {"B", 10.811},  // Boron
    {"C", 12.011},  // Carbon
    {"N", 14.007},  // Nitrogen
    {"O", 15.999},  // Oxygen
    {"F", 18.998},  // Fluorine
    {"Ne", 20.180}, // Neon
    {"Na", 22.990}, // Sodium
    {"Mg", 24.305}, // Magnesium
    {"Al", 26.982}, // Aluminum
    {"Si", 28.086}, // Silicon
    {"P", 30.974},  // Phosphorus
    {"S", 32.065},  // Sulfur
    {"Cl", 35.453}, // Chlorine
    {"K", 39.098},  // Potassium
    {"Ca", 40.078}, // Calcium
    {"Fe", 55.845}, // Iron
    {"Cu", 63.546}, // Copper
    {"Zn", 65.380}, // Zinc
    {"Br", 79.904}, // Bromine
    {"I", 126.904}, // Iodine
};

model::Structure PDBParser::parse_file(const std::string& filename) {
    model::Structure structure;
    if (!parse_to_structure(filename, structure)) {
        throw std::runtime_error("Failed to parse PDB file: " + filename);
    }
    return structure;
}

model::Structure PDBParser::parse_string(const std::string& pdbStr) {
    model::Structure structure;
    if (!parse_string_to_structure(pdbStr, structure)) {
        throw std::runtime_error("Failed to parse PDB string");
    }
    return structure;
}

bool PDBParser::parse_to_structure(const std::string& filename, model::Structure& structure) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Could not open PDB file: " << filename << std::endl;
        return false;
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return parse_string_to_structure(buffer.str(), structure);
}

bool PDBParser::parse_string_to_structure(const std::string& pdbStr, model::Structure& structure) {
    std::istringstream iss(pdbStr);
    std::string line;
    
    // Current residue tracking
    std::shared_ptr<model::Residue> currentResidue;
    
    // Clear existing data
    structure.clear();
    
    while (std::getline(iss, line)) {
        if (line.length() < 6) continue;
        
        RecordType recordType = getRecordType(line);
        bool success = true;
        
        switch (recordType) {
            case RecordType::ATOM:
            case RecordType::HETATM:
                success = parseAtomRecord(line, recordType, structure, currentResidue);
                break;
            case RecordType::TER:
                success = parseTerRecord(line, currentResidue, structure);
                break;
            case RecordType::HELIX:
                success = parseHelixRecord(line, structure);
                break;
            case RecordType::SHEET:
                success = parseSheetRecord(line, structure);
                break;
            case RecordType::SSBOND:
                success = parseSSBondRecord(line, structure);
                break;
            case RecordType::CRYST1:
                success = parseCryst1Record(line, structure);
                break;
            default:
                continue;
        }
        
        if (!success) {
            return false;
        }
    }
    
    // Calculate center of mass for the last residue if not already done
    if (currentResidue) {
        currentResidue->calculate_center_of_mass();
    }
    
    return true;
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

bool PDBParser::parseAtomRecord(const std::string& line, RecordType type,
                              model::Structure& structure,
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

        // Create atom with stripped values (using backward-compatible alias)
        auto atom = std::make_shared<model::Atom>();
        atom->set_bynu(serialNum);
        atom->set_type(atomName);
        atom->set_altloc(altLoc.empty() ? ' ' : altLoc[0]);
        atom->set_resname(resName);
        atom->set_chain(chainId.empty() ? ' ' : chainId[0]);
        atom->set_ires(resSeq);
        atom->set_inscode(iCode.empty() ? ' ' : iCode[0]);
        atom->set_coor(x, y, z);
        atom->set_occupancy(occupancy);
        atom->set_tempfactor(tempFactor);
        atom->set_segid(segId);
        atom->set_element(element);
        atom->set_charge_string(charge);
        atom->set_hetatm(type == RecordType::HETATM);

        // Set default mass based on element
        if (element.empty()) {
            // If element is not specified, guess from atom name
            std::string guessedElement = atomName;
            if (!guessedElement.empty()) {
                // Remove numbers and special characters
                guessedElement.erase(
                    std::remove_if(guessedElement.begin(), guessedElement.end(),
                                 [](char c) { return !std::isalpha(c); }),
                    guessedElement.end());
                if (!guessedElement.empty()) {
                    // Take first character and capitalize it
                    guessedElement = std::toupper(guessedElement[0]);
                    if (guessedElement.length() > 1) {
                        guessedElement += std::tolower(guessedElement[1]);
                    }
                }
            }
            element = guessedElement;
        }
        
        // Set mass based on element using the mass table
        double mass = 12.0; // Default to carbon mass if element not found
        auto it = ELEMENT_MASSES.find(element);
        if (it != ELEMENT_MASSES.end()) {
            mass = it->second;
        }
        atom->set_mass_charge(mass, 0.0); // Set mass and default charge to 0

        // Add atom to structure
        structure.add_atom(atom);

        // Handle residue
        if (!currentResidue || 
            currentResidue->get_chain() != chainId[0] ||
            currentResidue->get_ires() != resSeq ||
            currentResidue->get_inscode() != iCode[0]) {
            
            if (currentResidue) {
                currentResidue->calculate_center_of_mass();
            }
            
            // Create new residue
            currentResidue = std::make_shared<model::Residue>(
                resName, resSeq, segId, 0, chainId[0], iCode[0]);
            currentResidue->set_hetatm(type == RecordType::HETATM);
            structure.add_residue(currentResidue);
        }
        
        currentResidue->add_atom(atom);
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing ATOM/HETATM record: " << e.what() << std::endl;
        return false;
    }
}

bool PDBParser::parseTerRecord([[maybe_unused]] const std::string& line,
                             std::shared_ptr<model::Residue>& currentResidue,
                             model::Structure& structure) {
    try {
        if (currentResidue) {
            currentResidue->calculate_center_of_mass();
        }

        model::Structure::TerminalInfo terminal{
            currentResidue ? currentResidue->get_chain() : ' ',
            currentResidue ? currentResidue->get_ires() : 0,
            currentResidue ? currentResidue->get_inscode() : ' ',
            currentResidue ? currentResidue->get_resname() : ""
        };
        
        structure.add_terminal(terminal);
        currentResidue = nullptr;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing TER record: " << e.what() << std::endl;
        return false;
    }
}

bool PDBParser::parseHelixRecord(const std::string& line, model::Structure& structure) {
    try {
        model::Structure::SecondaryStructure helix;
        helix.id = line.substr(11, 3);
        helix.initResName = line.substr(15, 3);
        helix.initChainId = line[19];
        helix.initSeqNum = std::stoi(line.substr(21, 4));
        helix.initICode = line[25];
        helix.endResName = line.substr(27, 3);
        helix.endChainId = line[31];
        helix.endSeqNum = std::stoi(line.substr(33, 4));
        helix.endICode = line[37];
        helix.structureClass = std::stoi(line.substr(38, 2));
        
        std::string chainId(1, helix.initChainId);
        structure.add_helix(chainId, helix);
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing HELIX record: " << e.what() << std::endl;
        return false;
    }
}

bool PDBParser::parseSheetRecord(const std::string& line, model::Structure& structure) {
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
        
        structure.add_sheet(initChainId, sheetInfo);

        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing SHEET record: " << e.what() << std::endl;
        return false;
    }
}

bool PDBParser::parseSSBondRecord(const std::string& line, model::Structure& structure) {
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
        
        structure.add_ssbond(bondInfo);

        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing SSBOND record: " << e.what() << std::endl;
        return false;
    }
}

bool PDBParser::parseCryst1Record(const std::string& line, model::Structure& structure) {
    try {
        // Parse unit cell parameters according to PDB format
        double a = std::stod(line.substr(6, 9));
        double b = std::stod(line.substr(15, 9));
        double c = std::stod(line.substr(24, 9));
        double alpha = std::stod(line.substr(33, 7));
        double beta = std::stod(line.substr(40, 7));
        double gamma = std::stod(line.substr(47, 7));
        
        structure.set_box_dimensions(std::vector<double>{a, b, c, alpha, beta, gamma});
        
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing CRYST1 record: " << e.what() << std::endl;
        return false;
    }
}

} // namespace io
} // namespace pygcmc 