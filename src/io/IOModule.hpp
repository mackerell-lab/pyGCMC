#ifndef PYGCMC_IO_IOMODULE_HPP
#define PYGCMC_IO_IOMODULE_HPP

// Structure parsers
#include "structure/pdbParser.hpp"

// Topology parsers  
#include "topology/psfParser.hpp"
#include "topology/topParser.hpp"

// Force field parsers
#include "forcefield/prmParser.hpp"

// Parameter parsers
#include "parameters/inpParser.hpp"

namespace pygcmc {
namespace io {
    // Maintain backward compatibility - all original class names remain accessible
    using PDBParser = PDBParser;  // structure::PDBParser
    using PSFParser = PSFParser;  // topology::PSFParser  
    using TOPParser = TOPParser;  // topology::TOPParser
    using PRMParser = PRMParser;  // forcefield::PRMParser
    using INPParser = INPParser;  // parameters::INPParser
}
}

#endif // PYGCMC_IO_IOMODULE_HPP