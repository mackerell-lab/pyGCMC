#pragma once

#include "../../../../model/ModelModule.hpp"
#include "../../movement/reservoir/fragment_reservoir.hpp"
#include <memory>
#include <string>
#include <filesystem>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {
namespace setup {

/**
 * @brief Builds complete simulation input from various file sources
 * 
 * This class orchestrates the loading of INP, PDB, TOP, PAR, and ITP files
 * to create a fully initialized MCState with real molecular structure,
 * topology, and force field parameters.
 */
class SimulationInputBuilder {
public:
    /**
     * @brief Configuration for input building
     */
    struct Config {
        std::string inpFile;           // INP parameter file
        bool loadStructure = true;     // Load PDB structure if available
        bool loadTopology = true;      // Load TOP topology if available
        bool loadParameters = true;    // Load PAR force field parameters
        bool loadFragments = true;     // Load ITP fragment templates
        bool verbose = false;          // Verbose logging
    };

    /**
     * @brief Result of building simulation input
     */
    struct Result {
        std::shared_ptr<model::param::Param> parameters;        // Parsed INP parameters
        std::shared_ptr<model::Molecular> molecular;            // Combined molecular system
        std::shared_ptr<model::ForceField> forceField;         // Force field parameters
        std::shared_ptr<model::montecarlo::MCState> mcState;   // Initialized MC state
        std::map<std::string, platform::cpu::movement::FragmentTemplate> fragmentTemplates; // Fragment templates
        bool structureLoaded = false;                          // Whether PDB was loaded
        bool topologyLoaded = false;                           // Whether TOP was loaded
        bool parametersLoaded = false;                         // Whether PAR was loaded
    };

    SimulationInputBuilder(const Config& config);
    ~SimulationInputBuilder() = default;

    /**
     * @brief Build complete simulation input
     * 
     * This method:
     * 1. Parses INP file to get all file paths and parameters
     * 2. Loads PDB structure if specified
     * 3. Loads TOP topology if specified
     * 4. Loads PAR force field parameters if specified
     * 5. Combines molecular data using MolecularCombiner
     * 6. Initializes MCState using MCInitializer
     * 
     * @return Result containing all loaded data
     * @throw std::runtime_error if required files cannot be loaded
     */
    Result build();

    /**
     * @brief Load fragment templates from ITP files
     * 
     * @param fragItpFiles List of ITP file paths from INP
     * @param parameters Parsed parameters containing fragment info
     * @return Map of fragment name to template
     */
    std::map<std::string, platform::cpu::movement::FragmentTemplate> loadFragmentTemplates(
        const std::vector<std::string>& fragItpFiles,
        const std::shared_ptr<model::param::Param>& parameters);

private:
    Config config_;

    /**
     * @brief Parse INP file
     */
    std::shared_ptr<model::param::Param> parseINP(const std::string& filename);

    /**
     * @brief Load PDB structure
     */
    std::shared_ptr<model::Structure> loadPDB(const std::string& filename);

    /**
     * @brief Load TOP topology
     */
    std::shared_ptr<model::Topology> loadTopology(const std::string& filename);

    /**
     * @brief Load PAR force field parameters
     */
    std::shared_ptr<model::ForceField> loadParameters(const std::vector<std::string>& parFiles);

    /**
     * @brief Combine molecular data
     */
    std::shared_ptr<model::Molecular> combineMolecular(
        const std::shared_ptr<model::Structure>& structure,
        const std::shared_ptr<model::Topology>& topology,
        const std::shared_ptr<model::ForceField>& forceField);

    /**
     * @brief Initialize MC state from molecular data
     */
    std::shared_ptr<model::MCState> initializeMCState(
        const std::shared_ptr<model::Molecular>& molecular,
        const std::shared_ptr<model::ForceField>& forceField);

    /**
     * @brief Resolve file path relative to base directory
     */
    std::string resolveFilePath(const std::string& path, 
                                const std::filesystem::path& baseDir) const;
    
    void log(const std::string& message) const;
};

} // namespace setup
} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc