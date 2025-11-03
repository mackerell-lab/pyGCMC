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
#include <filesystem>

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
    
    // Determine base directory for resolving relative paths
    std::filesystem::path baseDir = std::filesystem::path(config_.inpFile).parent_path();
    if (baseDir.empty()) {
        baseDir = std::filesystem::current_path();
    }
    
    const auto& fileInfo = result.parameters->get_file_info();
    
    // Step 2: Load PDB structure if available
    std::shared_ptr<Structure> structure;
    if (config_.loadStructure && !fileInfo.input_pdb_file.empty()) {
        std::string pdbPath = resolveFilePath(fileInfo.input_pdb_file, baseDir);
        log("Loading PDB structure: " + pdbPath);
        try {
            structure = loadPDB(pdbPath);
            if (structure) {
                result.structureLoaded = true;
                log("Loaded " + std::to_string(structure->get_atoms().size()) + " atoms from PDB");
            }
        } catch (const std::exception& e) {
            log("ERROR: Failed to load PDB: " + std::string(e.what()));
            throw;
        }
    }
    
    // Step 3: Load TOP topology if available
    std::shared_ptr<Topology> topology;
    if (config_.loadTopology && !fileInfo.topology_file.empty()) {
        std::string topPath = resolveFilePath(fileInfo.topology_file, baseDir);
        log("Loading topology: " + topPath);
        try {
            topology = loadTopology(topPath);
            if (topology) {
                result.topologyLoaded = true;
                log("Loaded topology with " + std::to_string(topology->get_num_atoms()) + " atoms");
            }
        } catch (const std::exception& e) {
            log("ERROR: Failed to load topology: " + std::string(e.what()));
            throw;
        }
    }
    
    // Step 4: Load PAR force field parameters if available
    if (config_.loadParameters && !fileInfo.par_files.empty()) {
        log("Loading force field parameters");
        try {
            std::vector<std::string> resolvedParFiles;
            for (const auto& parFile : fileInfo.par_files) {
                resolvedParFiles.push_back(resolveFilePath(parFile, baseDir));
            }
            result.forceField = loadParameters(resolvedParFiles);
            if (result.forceField) {
                result.parametersLoaded = true;
                log("Loaded force field parameters");
            }
        } catch (const std::exception& e) {
            log("Warning: Failed to load parameters: " + std::string(e.what()));
        }
    }
    
    // Step 4b: Load fragment templates if available
    if (!fileInfo.fragment_top_files.empty()) {
        log("Loading fragment templates");
        std::vector<std::string> resolvedFragFiles;
        for (const auto& fragFile : fileInfo.fragment_top_files) {
            resolvedFragFiles.push_back(resolveFilePath(fragFile, baseDir));
        }
        result.fragmentTemplates = loadFragmentTemplates(resolvedFragFiles, 
                                                         result.parameters);
        log("Loaded " + std::to_string(result.fragmentTemplates.size()) + " fragment templates");
    }
    
    // Step 5: Combine molecular data if we have BOTH structure and topology with actual data
    if (result.structureLoaded && result.topologyLoaded) {
        log("Combining molecular data");
        // Use already loaded structure instead of re-loading
        
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
        result.molecular->atoms.size() > 0 && result.forceField) {
        // Use MCInitializer to populate from real molecular data
        system::montecarlo::MCInitializer initializer;
        initializer.initializeFromMolecular(*result.mcState, result.molecular);
        
        initializer.initializeForceField(*result.mcState, *result.forceField);
        
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

        // Set cutoff from space info
        result.mcState->info.cutoff = spaceInfo.cutoff;

        // Also sync cutoff to energy parameters if parameters exist
        if (result.parameters) {
            auto& energyInfo = const_cast<model::param::EnergyInfo&>(result.parameters->get_energy_info());
            energyInfo.fragment_cutoff = spaceInfo.cutoff;
            energyInfo.protein_cutoff = spaceInfo.cutoff;
            energyInfo.fragment_cutoff_squared = spaceInfo.cutoff * spaceInfo.cutoff;
            energyInfo.protein_cutoff_squared = spaceInfo.cutoff * spaceInfo.cutoff;
        }

        log("MC state initialized with INP parameters only");
    }
    
    return result;
}

std::map<std::string, platform::cpu::movement::FragmentTemplate> SimulationInputBuilder::loadFragmentTemplates(
    const std::vector<std::string>& fragItpFiles,
    const std::shared_ptr<model::param::Param>& parameters) {
    
    std::map<std::string, platform::cpu::movement::FragmentTemplate> templates;
    
    if (!parameters) {
        throw std::runtime_error("Parameters not available for fragment template loading");
    }
    
    // Get fragment info from parameters for matching
    const auto& fragNames = parameters->get_file_info().fragment_names;
    const auto& fragConcs = parameters->get_fragment_info().conc_list;
    const auto& fragMuexs = parameters->get_fragment_info().muex_list;
    
    for (size_t i = 0; i < fragItpFiles.size(); ++i) {
        const auto& itpFile = fragItpFiles[i];
        log("Loading fragment template from: " + itpFile);
        
        // Check if file exists
        std::ifstream file(itpFile);
        if (!file.good()) {
            log("ERROR: Fragment template not found: " + itpFile);
            throw std::runtime_error("Fragment template not found: " + itpFile);
        }
        file.close();
        
        // Extract fragment name from filename
        std::filesystem::path itpPath(itpFile);
        std::string fragName = itpPath.stem().string();
        
        // Use FragmentLibrary to load the ITP file
        io::topology::FragmentLibrary fragLib;
        if (!fragLib.loadFromITP(itpFile, fragName, templates.size())) {
            throw std::runtime_error("Failed to parse fragment template: " + itpFile);
        }
        
        // Get the fragment data
        auto fragmentData = fragLib.get(fragName);
        if (!fragmentData) {
            throw std::runtime_error("Failed to retrieve fragment data for: " + fragName);
        }
        
        // Build complete FragmentTemplate
        platform::cpu::movement::FragmentTemplate tmpl;
        tmpl.name = fragmentData->name;
        tmpl.typeId = fragmentData->typeId;
        tmpl.atoms = fragmentData->atoms;
        tmpl.molecularWeight = fragmentData->molecularWeight;
        tmpl.radius = fragmentData->radius;
        
        // Find matching fragment in parameters to get concentration and chemical potential
        auto nameIt = std::find(fragNames.begin(), fragNames.end(), fragName);
        if (nameIt != fragNames.end()) {
            size_t idx = std::distance(fragNames.begin(), nameIt);
            if (idx < fragConcs.size()) {
                tmpl.concentration = fragConcs[idx];
            }
            if (idx < fragMuexs.size()) {
                tmpl.chemicalPotential = fragMuexs[idx];
                // Calculate activity: z = exp(β*μ)
                double beta = parameters->get_mc_info().beta;
                tmpl.activity = std::exp(beta * tmpl.chemicalPotential);
            }
        }
        
        templates[tmpl.name] = tmpl;
        log("Loaded fragment " + tmpl.name + " with " + 
            std::to_string(tmpl.atoms.size()) + " atoms, " +
            "conc=" + std::to_string(tmpl.concentration) + " M, " +
            "μ=" + std::to_string(tmpl.chemicalPotential) + " kJ/mol");
    }
    
    if (templates.empty() && !fragItpFiles.empty()) {
        throw std::runtime_error("Failed to load any fragment templates");
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

std::string SimulationInputBuilder::resolveFilePath(const std::string& path, 
                                                    const std::filesystem::path& baseDir) const {
    std::filesystem::path filePath(path);
    
    // If path is already absolute, return as-is
    if (filePath.is_absolute()) {
        return path;
    }
    
    // Otherwise resolve relative to base directory
    std::filesystem::path resolved = baseDir / filePath;
    return resolved.string();
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
