#ifndef PYGCMC_IO_IOMODULE_HPP
#define PYGCMC_IO_IOMODULE_HPP

// Structure parsers
#include "structure/PdbParserMain.hpp"

// Topology parsers  
#include "topology/psfParser.hpp"
#include "topology/topParser.hpp"

// Force field parsers
#include "forcefield/prmParser.hpp"

// Parameter parsers
#include "parameters/InpParserMain.hpp"

namespace pygcmc {
namespace io {
    // Maintain backward compatibility - all original class names remain accessible
    using PDBParser = structure::PdbParserMain;  // structure::PdbParserMain
    using PSFParser = PSFParser;  // topology::PSFParser  
    using TOPParser = TOPParser;  // topology::TOPParser
    using PRMParser = PRMParser;  // forcefield::PRMParser
    using INPParser = parameters::InpParserMain;  // parameters::InpParserMain
}
}

#endif // PYGCMC_IO_IOMODULE_HPP