#ifndef PYGCMC_IO_IOMODULE_HPP
#define PYGCMC_IO_IOMODULE_HPP

// Structure parsers
#include "structure/PdbParserMain.hpp"

// Topology parsers  
#include "topology/psfParser.hpp"
#include "topology/topParser.hpp"

// Force field parsers
#include "forcefield/PrmParserMain.hpp"

// Parameter parsers
#include "parameters/InpParserMain.hpp"

namespace pygcmc {
namespace io {
    // Maintain backward compatibility - all original class names remain accessible
    using PDBParser = structure::PdbParserMain;  // Map structure::PdbParserMain to PDBParser
    using INPParser = parameters::InpParserMain;  // Map parameters::InpParserMain to INPParser
    
    // Note: PSFParser, TOPParser, and PRMParser are already defined directly in pygcmc::io namespace
    // so no using declarations are needed for them
}
}

#endif // PYGCMC_IO_IOMODULE_HPP