// src/io/structure/PdbParserRecords.cpp

#include "PdbParserRecords.hpp"
#include "PdbParserStructureRecords.hpp"
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <set>

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
        
        bool needNewResidue = false;
        
        // Standard residue identification criteria
        if (!currentResidue || currentResidue->get_resname() != resName || 
            currentResidue->get_chain() != chainChar || 
            currentResidue->get_ires() != resSeq || 
            currentResidue->get_inscode() != iCodeChar) {
            needNewResidue = true;
        }
        // Additional continuity check for small molecules - non-consecutive atoms = different molecules
        else if (PdbParserStructureRecords::isSmallMolecule(resName) && currentResidue->atom_count() > 0) {
            // Get the last atom's serial number from current residue
            const auto& atoms = currentResidue->get_atoms();
            if (!atoms.empty()) {
                int lastSerialNum = atoms.back()->get_bynu();
                // If current atom serial number is not consecutive, it's a different molecule
                if (serialNum != lastSerialNum + 1) {
                    needNewResidue = true;
                }
            }
        }
        
        if (needNewResidue) {
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

} // namespace structure
} // namespace io
} // namespace pygcmc