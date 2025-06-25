// src/io/structure/PdbParserStructures.cpp

#include "PdbParserStructures.hpp"

namespace pygcmc {
namespace io {
namespace structure {

// Define the element mass table
const std::map<std::string, double> PdbParserStructures::ELEMENT_MASSES = {
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

PdbParserStructures::RecordType PdbParserStructures::getRecordType(const std::string& line) {
    if (line.length() < 6) return RecordType::UNKNOWN;
    
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

const std::map<std::string, double>& PdbParserStructures::getElementMasses() {
    return ELEMENT_MASSES;
}

} // namespace structure
} // namespace io
} // namespace pygcmc