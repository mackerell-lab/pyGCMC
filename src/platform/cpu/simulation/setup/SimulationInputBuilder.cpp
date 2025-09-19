#include "SimulationInputBuilder.hpp"
#include "../../../../io/parameters/InpParserGCMC.hpp"
#include "../../../../io/structure/PdbParserMain.hpp"
#include "../../../../io/topology/topParserMain.hpp"
#include "../../../../io/forcefield/PrmParserMain.hpp"
#include "../../../../io/topology/FragmentLibrary.hpp"
#include "../../../../system/molecular/MolecularCombiner.hpp"
#include "../../../../system/montecarlo/MCInitializer.hpp"
#include "../../../../system/log/LogMain.hpp"
#include <fstream>
#include <stdexcept>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {
namespace setup {

using namespace pygcmc::model;
using namespace pygcmc::io;

SimulationInputBuilder::SimulationInputBuilder(const Config& config)
    : config_(config) {
}

SimulationInputBuilder::Result SimulationInputBuilder::build() {
    Result result;
    
    // Step 1: Parse INP file
    log("Parsing INP file: " + config_.inpFile);
    result.parameters = parseINP(config_.inpFile);
    if (!result.parameters) {
        throw std::runtime_error("Failed to parse INP file");
    }
    
    const auto& fileInfo = result.parameters->get_file_info();
    
    // Step 2: Load PDB structure if available
    if (config_.loadStructure && !fileInfo.input_pdb_file.empty()) {
        log("Loading PDB structure: " + fileInfo.input_pdb_file);
        try {
            auto structure = loadPDB(fileInfo.input_pdb_file);
            if (structure) {
                result.structureLoaded = true;
                log("Loaded " + std::to_string(structure->get_atoms().size()) + " atoms from PDB");
            }
        } catch (const std::exception& e) {
            log("Warning: Failed to load PDB: " + std::string(e.what()));
        }
    }
    
    // Step 3: Load TOP topology if available
    std::shared_ptr<Topology> topology;
    if (config_.loadTopology && !fileInfo.topology_file.empty()) {
        log("Loading topology: " + fileInfo.topology_file);
        try {
            topology = loadTopology(fileInfo.topology_file);
            if (topology) {
                result.topologyLoaded = true;
                log("Loaded topology with " + std::to_string(topology->get_num_atoms()) + " atoms");
            }
        } catch (const std::exception& e) {
            log("Warning: Failed to load topology: " + std::string(e.what()));
        }
    }
    
    // Step 4: Load PAR force field parameters if available
    if (config_.loadParameters && !fileInfo.par_files.empty()) {
        log("Loading force field parameters");
        try {
            result.forceField = loadParameters(fileInfo.par_files);
            if (result.forceField) {
                result.parametersLoaded = true;
                log("Loaded force field parameters");
            }
        } catch (const std::exception& e) {
            log("Warning: Failed to load parameters: " + std::string(e.what()));
        }
    }
    
    // Step 5: Combine molecular data if we have BOTH structure and topology with actual data
    if (result.structureLoaded && result.topologyLoaded) {
        log("Combining molecular data");
        auto structure = loadPDB(fileInfo.input_pdb_file);
        
        // Only combine if we have actual atoms
        if (structure && structure->get_atoms().size() > 0 && 
            topology && topology->get_num_atoms() > 0) {
            result.molecular = combineMolecular(structure, topology, result.forceField);
        } else {
            log("Warning: Structure or topology is empty, skipping molecular combination");
            result.structureLoaded = false;  // Mark as not loaded if empty
            result.topologyLoaded = false;
        }
    }
    
    // Step 6: Initialize MC state
    log("Initializing MC state");
    result.mcState = std::make_shared<montecarlo::MCState>();
    
    // Only use MCInitializer if we have complete molecular data with actual atoms
    if (result.molecular && result.structureLoaded && result.topologyLoaded &&
        result.molecular->atoms.size() > 0) {
        // Use MCInitializer to populate from real molecular data
        system::montecarlo::MCInitializer initializer;
        initializer.initializeFromMolecular(*result.mcState, result.molecular);
        
        if (result.forceField) {
            initializer.initializeForceField(*result.mcState, *result.forceField);
        }
        
        log("MC state initialized with real molecular data");
    } else {
        // Fall back to empty state with box dimensions from INP
        const auto& spaceInfo = result.parameters->get_space_info();
        result.mcState->info.box[0] = spaceInfo.box_size[0];
        result.mcState->info.box[1] = spaceInfo.box_size[1];
        result.mcState->info.box[2] = spaceInfo.box_size[2];
        
        // Set temperature
        const auto& mcInfo = result.parameters->get_mc_info();
        result.mcState->info.beta = mcInfo.beta;
        
        log("MC state initialized with INP parameters only");
    }
    
    return result;
}

std::map<std::string, platform::cpu::movement::FragmentTemplate> SimulationInputBuilder::loadFragmentTemplates(
    const std::vector<std::string>& fragItpFiles) {
    
    std::map<std::string, platform::cpu::movement::FragmentTemplate> templates;
    
    for (const auto& itpFile : fragItpFiles) {
        try {
            log("Loading fragment template from: " + itpFile);
            
            // Parse ITP file (implementation depends on ITPParser availability)
            // For now, create placeholder
            // TODO: Implement proper ITP parsing
            
            // Extract fragment name from filename
            size_t lastSlash = itpFile.find_last_of("/\\");
            size_t lastDot = itpFile.find_last_of(".");
            std::string fragName = itpFile.substr(
                lastSlash != std::string::npos ? lastSlash + 1 : 0,
                lastDot - (lastSlash != std::string::npos ? lastSlash + 1 : 0)
            );
            
            platform::cpu::movement::FragmentTemplate tmpl;
            tmpl.name = fragName;
            // TODO: Populate from ITP
            
            templates[fragName] = tmpl;
            
        } catch (const std::exception& e) {
            log("Warning: Failed to load fragment template from " + itpFile + ": " + e.what());
        }
    }
    
    return templates;
}

std::shared_ptr<param::Param> SimulationInputBuilder::parseINP(const std::string& filename) {
    auto params = std::make_shared<param::Param>();
    parameters::InpParserGCMC::parse_to_param(filename, *params);
    parameters::InpParserGCMC::enhance_param(*params);
    return params;
}

std::shared_ptr<model::Structure> SimulationInputBuilder::loadPDB(const std::string& filename) {
    // Check if file exists
    std::ifstream file(filename);
    if (!file.good()) {
        throw std::runtime_error("PDB file not found: " + filename);
    }
    file.close();
    
    auto structure = std::make_shared<model::Structure>();
    io::structure::PdbParserMain::parse_to_structure(filename, *structure);
    return structure;
}

std::shared_ptr<model::Topology> SimulationInputBuilder::loadTopology(const std::string& filename) {
    // Check if file exists
    std::ifstream file(filename);
    if (!file.good()) {
        throw std::runtime_error("Topology file not found: " + filename);
    }
    file.close();
    
    auto topology = std::make_shared<model::Topology>();
    io::TOPParser parser;
    parser.parse_to_topology(filename, *topology);
    return topology;
}

std::shared_ptr<model::ForceField> SimulationInputBuilder::loadParameters(
    const std::vector<std::string>& parFiles) {
    
    auto forceField = std::make_shared<model::ForceField>();
    
    for (const auto& parFile : parFiles) {
        try {
            log("Loading parameters from: " + parFile);
            io::PRMParser::parse_file_to_forcefield(parFile, *forceField);
        } catch (const std::exception& e) {
            log("Warning: Failed to load " + parFile + ": " + e.what());
        }
    }
    
    return forceField;
}

std::shared_ptr<model::Molecular> SimulationInputBuilder::combineMolecular(
    const std::shared_ptr<model::Structure>& structure,
    const std::shared_ptr<model::Topology>& topology,
    const std::shared_ptr<model::ForceField>& /*forceField*/) {
    
    system::molecular::MolecularCombiner combiner;
    
    // Create molecular system
    auto molecular = std::make_shared<model::Molecular>();
    
    if (structure && topology) {
        // Combine structure and topology
        molecular = combiner.combine(structure, topology);
    } else if (structure) {
        // Structure only - copy atoms and residues
        molecular->atoms = structure->get_atoms();
        molecular->residues = structure->get_residues();
    } else if (topology) {
        // Topology only - copy available topology data
        // Note: Topology class doesn't expose atoms/residues directly
        molecular->bonds = topology->get_bonds();
        molecular->angles = topology->get_angles();
        molecular->dihedrals = topology->get_dihedrals();
    }
    
    // Note: ForceField is separate and not stored in Molecular
    
    return molecular;
}

std::shared_ptr<model::MCState> SimulationInputBuilder::initializeMCState(
    const std::shared_ptr<model::Molecular>& molecular,
    const std::shared_ptr<model::ForceField>& forceField) {
    
    auto mcState = std::make_shared<model::MCState>();
    
    if (molecular) {
        system::montecarlo::MCInitializer initializer;
        initializer.initializeFromMolecular(*mcState, molecular);
        
        if (forceField) {
            initializer.initializeForceField(*mcState, *forceField);
        }
    }
    
    return mcState;
}

void SimulationInputBuilder::log(const std::string& message) const {
    if (config_.verbose) {
        system::log::LogMain::info("[SimulationInputBuilder] ", message);
    }
}

} // namespace setup
} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc