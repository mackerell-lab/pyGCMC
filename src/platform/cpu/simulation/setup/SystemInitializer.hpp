#ifndef PYGCMC_PLATFORM_CPU_SIMULATION_SETUP_SYSTEM_INITIALIZER_HPP
#define PYGCMC_PLATFORM_CPU_SIMULATION_SETUP_SYSTEM_INITIALIZER_HPP

#include "../../movement/gcmc/GCMCEngine.hpp"
#include "../../movement/gcmc/GCMCAcceptance.hpp"
#include "../../movement/reservoir/MultiTypeReservoir.hpp"
#include "../../movement/reservoir/fragment_reservoir.hpp"
#include "../../../../model/param/ParamMain.hpp"
#include "../../../../model/montecarlo/MCMain.hpp"
#include <memory>
#include <vector>
#include <string>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {

/**
 * @brief System setup and initialization module
 *
 * This class handles all initialization tasks including:
 * - Parameter loading
 * - System setup
 * - Fragment configuration
 * - Engine initialization
 */
class SystemInitializer {
public:
    /**
     * @brief Fragment configuration
     */
    struct FragmentConfig {
        std::string name;
        int typeId;
        double concentration;
        double chemicalPotential;
        double activity;
        double probability;
        int maxCount;
        movement::FragmentTemplate template_;
    };

    // Constructor
    SystemInitializer();
    ~SystemInitializer();

    // Parameter loading
    bool loadParameters(const std::string& inputFile,
                       model::param::Param& params);

    // System setup
    bool setupSystem(const model::param::Param& params,
                    model::montecarlo::MCState& state);

    // Fragment setup
    bool setupFragments(const model::param::Param& params,
                       std::vector<FragmentConfig>& fragments,
                       movement::MultiTypeReservoir& reservoir);

    // Engine setup
    bool setupEngine(movement::gcmc::GCMCEngine& engine,
                    model::montecarlo::MCState* state,
                    movement::MultiTypeReservoir* reservoir,
                    const model::param::Param& params);

    // Acceptance calculator setup
    bool setupAcceptance(movement::gcmc::GCMCAcceptance& acceptance,
                        const model::param::Param& params);

    // Validation
    bool validateSetup(const model::param::Param& params,
                      const model::montecarlo::MCState& state);

private:
    // Helper methods
    bool loadTopology(const std::string& filename,
                     model::montecarlo::MCState& state);

    bool loadCoordinates(const std::string& filename,
                        model::montecarlo::MCState& state);

    bool loadForceField(const std::string& filename,
                       model::montecarlo::MCForceField& ff);

    bool loadFragmentTemplate(const std::string& filename,
                             movement::FragmentTemplate& tmpl);

    void calculateActivities(std::vector<FragmentConfig>& fragments,
                           double temperature);

    void calculateProbabilities(std::vector<FragmentConfig>& fragments);

    // Error handling
    void reportError(const std::string& message) const;
};

} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_SIMULATION_SETUP_SYSTEM_INITIALIZER_HPP
