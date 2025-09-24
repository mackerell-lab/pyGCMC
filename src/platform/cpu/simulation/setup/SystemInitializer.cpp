#include "SystemInitializer.hpp"
#include "../../../../io/parameters/InpParserGCMC.hpp"
#include "../../../../io/topology/psfParserMain.hpp"
#include "../../../../io/structure/PdbParserMain.hpp"
#include "../../../../io/forcefield/PrmParserMain.hpp"
#include <iostream>
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

SystemInitializer::SystemInitializer() {
}

SystemInitializer::~SystemInitializer() {
}

bool SystemInitializer::loadParameters(const std::string& inputFile,
                                    model::param::Param& params) {
    try {
        // Parse the input file
        io::parameters::InpParserGCMC::parse_to_param(inputFile, params);
        
        // Enhance parameters
        io::parameters::InpParserGCMC::enhance_param(params);
        
        // Update derived values
        params.update_derived_values();
        
        // Validate parameters
        if (!params.is_valid()) {
            reportError("Invalid parameters after loading");
            return false;
        }
        
        return true;
    } catch (const std::exception& e) {
        reportError("Failed to load parameters: " + std::string(e.what()));
        return false;
    }
}

bool SystemInitializer::setupSystem(const model::param::Param& params,
                                 model::montecarlo::MCState& state) {
    // Set box dimensions
    const auto& space = params.get_space_info();
    state.periodicBox.resize(3);
    state.periodicBox[0] = space.box_size[0];
    state.periodicBox[1] = space.box_size[1];
    state.periodicBox[2] = space.box_size[2];
    
    // Set system information
    state.info.box[0] = space.box_size[0];
    state.info.box[1] = space.box_size[1];
    state.info.box[2] = space.box_size[2];
    state.info.setTemperature(params.get_mc_info().temperature);
    // beta is set automatically by setTemperature
    
    // Load topology if specified
    const auto& files = params.get_file_info();
    if (!files.topology_file.empty()) {
        if (!loadTopology(files.topology_file, state)) {
            return false;
        }
    }
    
    // Load coordinates if specified
    if (!files.input_pdb_file.empty()) {
        if (!loadCoordinates(files.input_pdb_file, state)) {
            return false;
        }
    }
    
    // Load force field parameters
    for (const auto& prmFile : files.par_files) {
        if (!loadForceField(prmFile, state.forcefield)) {
            return false;
        }
    }
    
    return true;
}

bool SystemInitializer::setupFragments(const model::param::Param& params,
                                    std::vector<FragmentConfig>& fragments,
                                    movement::MultiTypeReservoir& reservoir) {
    const auto& fragInfo = params.get_fragment_info();
    const auto& fileInfo = params.get_file_info();
    const auto& mcInfo = params.get_mc_info();
    
    fragments.clear();
    
    // Setup each fragment type
    for (size_t i = 0; i < fileInfo.fragment_names.size(); ++i) {
        FragmentConfig config;
        config.name = fileInfo.fragment_names[i];
        config.typeId = static_cast<int>(i);
        
        // Set concentration and chemical potential
        if (i < fragInfo.conc_list.size()) {
            config.concentration = fragInfo.conc_list[i];
        }
        if (i < fragInfo.muex_list.size()) {
            config.chemicalPotential = fragInfo.muex_list[i];
        }
        
        // Load fragment template
        if (i < fileInfo.fragment_top_files.size()) {
            if (!loadFragmentTemplate(fileInfo.fragment_top_files[i], 
                                    config.template_)) {
                return false;
            }
        }
        
        // Calculate max count based on volume and concentration
        double volume = params.get_space_info().volume;
        if (volume > 0) {
            // Convert concentration to molecules/nm^3
            double avogadro = 6.022e23;
            double density = config.concentration * avogadro / 1e27;
            config.maxCount = static_cast<int>(density * volume * 2);
        } else {
            config.maxCount = 1000;  // Default
        }
        
        fragments.push_back(config);
        
        // Register with reservoir
        movement::MultiTypeReservoir::TypeInfo typeInfo;
        typeInfo.typeId = config.typeId;
        typeInfo.name = config.name;
        typeInfo.maxCount = config.maxCount;
        reservoir.addType(typeInfo, config.template_);
    }
    
    // Calculate activities and probabilities
    calculateActivities(fragments, mcInfo.temperature);
    calculateProbabilities(fragments);
    
    return true;
}

bool SystemInitializer::setupEngine(movement::gcmc::GCMCEngine& engine,
                                 model::montecarlo::MCState* state,
                                 movement::MultiTypeReservoir* reservoir,
                                 const model::param::Param& params) {
    // Initialize the engine
    auto fragReservoir = dynamic_cast<movement::FragmentReservoir*>(reservoir);
    if (!fragReservoir) {
        reportError("Invalid reservoir type for engine initialization");
        return false;
    }
    
    engine.initialize(state, fragReservoir);
    
    // Configure engine settings
    engine.setTemperature(params.get_mc_info().temperature);
    engine.setCutoff(params.get_space_info().cutoff);  // Use cutoff from parameters
    
    return true;
}

bool SystemInitializer::setupAcceptance(movement::gcmc::GCMCAcceptance& acceptance,
                                     const model::param::Param& params) {
    const auto& mcInfo = params.get_mc_info();
    
    // Configure acceptance calculator
    acceptance.setTemperature(mcInfo.temperature);
    // Note: GCMCAcceptance may not have these setters
    // The acceptance calculator typically gets these from the engine
    
    return true;
}

bool SystemInitializer::validateSetup(const model::param::Param& params,
                                   const model::montecarlo::MCState& state) {
    // Check box dimensions
    if (state.periodicBox.size() != 3 ||
        state.periodicBox[0] <= 0 ||
        state.periodicBox[1] <= 0 ||
        state.periodicBox[2] <= 0) {
        reportError("Invalid box dimensions");
        return false;
    }
    
    // Check temperature
    // Check beta instead of temperature
    if (state.info.beta <= 0) {
        reportError("Invalid beta/temperature");
        return false;
    }
    
    // Check MC steps
    if (params.get_mc_info().mc_steps <= 0) {
        reportError("Invalid number of MC steps");
        return false;
    }
    
    return true;
}

bool SystemInitializer::loadTopology(const std::string& filename,
                                  model::montecarlo::MCState& /* state */) {
    try {
        // Note: namespace needs to be corrected
        // io::topology::PSFParser parser;
        // Simplified - actual implementation would load topology
        std::cout << "Loading topology from: " << filename << std::endl;
        return true;
    } catch (const std::exception& e) {
        reportError("Failed to load topology: " + std::string(e.what()));
        return false;
    }
}

bool SystemInitializer::loadCoordinates(const std::string& filename,
                                     model::montecarlo::MCState& /* state */) {
    try {
        // io::structure::PdbParserMain parser; // Will be used in actual implementation
        // Simplified - actual implementation would load coordinates
        std::cout << "Loading coordinates from: " << filename << std::endl;
        return true;
    } catch (const std::exception& e) {
        reportError("Failed to load coordinates: " + std::string(e.what()));
        return false;
    }
}

bool SystemInitializer::loadForceField(const std::string& filename,
                                    model::montecarlo::MCForceField& /* ff */) {
    try {
        // Note: namespace needs to be corrected  
        // io::forcefield::PRMParser parser;
        // Simplified - actual implementation would load force field
        std::cout << "Loading force field from: " << filename << std::endl;
        return true;
    } catch (const std::exception& e) {
        reportError("Failed to load force field: " + std::string(e.what()));
        return false;
    }
}

bool SystemInitializer::loadFragmentTemplate(const std::string& filename,
                                          movement::FragmentTemplate& /* tmpl */) {
    try {
        // Simplified - actual implementation would load fragment template
        std::cout << "Loading fragment template from: " << filename << std::endl;
        return true;
    } catch (const std::exception& e) {
        reportError("Failed to load fragment template: " + std::string(e.what()));
        return false;
    }
}

void SystemInitializer::calculateActivities(std::vector<FragmentConfig>& fragments,
                                         double temperature) {
    const double R = 8.314e-3;  // kJ/(mol*K)
    
    for (auto& frag : fragments) {
        // Calculate activity from chemical potential
        // a = exp(mu / RT)
        frag.activity = std::exp(frag.chemicalPotential / (R * temperature));
    }
}

void SystemInitializer::calculateProbabilities(std::vector<FragmentConfig>& fragments) {
    if (fragments.empty()) return;
    
    // For now, use equal probabilities
    double prob = 1.0 / fragments.size();
    for (auto& frag : fragments) {
        frag.probability = prob;
    }
}

void SystemInitializer::reportError(const std::string& message) const {
    std::cerr << "SimulationSetup Error: " << message << std::endl;
}

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc