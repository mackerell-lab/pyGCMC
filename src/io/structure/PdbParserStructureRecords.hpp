// src/io/structure/PdbParserStructureRecords.hpp

#pragma once

#include "PdbParserStructures.hpp"
#include "model/ModelModule.hpp"
#include <string>
#include <memory>

namespace pygcmc {
namespace io {
namespace structure {

class PdbParserStructureRecords {
public:
    // Parse secondary structure and crystal info records
    static bool parseHelixRecord(const std::string& line, model::Structure& structure);
    static bool parseSheetRecord(const std::string& line, model::Structure& structure);
    static bool parseSSBondRecord(const std::string& line, model::Structure& structure);
    static bool parseCryst1Record(const std::string& line, model::Structure& structure);
    
    // Helper function for small molecule identification
    static bool isSmallMolecule(const std::string& resName);
};

} // namespace structure
} // namespace io
} // namespace pygcmc