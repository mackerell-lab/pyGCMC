#pragma once

/**
 * @file GCMCController.hpp
 * @brief GCMC-specific control logic
 */

#include "../../../../model/montecarlo/MCMain.hpp"
#include "../../movement/gcmc/GCMCEngine.hpp"
#include "../../movement/gcmc/GCMCAcceptance.hpp"
#include "../../movement/gcmc/GCMCStatistics.hpp"
#include "../../movement/reservoir/MultiTypeReservoir.hpp"
#include <memory>
#include <map>
#include <string>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace simulation {
namespace gcmc {

/**
 * @brief Controls GCMC-specific operations and move selection
 */
class GCMCController {
public:
    /**
     * @brief Configuration for GCMC control
     */
    struct Config {
        // Move probabilities (should sum to 1.0)
        double insertionProb = 0.25;
        double deletionProb = 0.25;
        double translationProb = 0.25;
        double rotationProb = 0.25;
        
        // Advanced options
        bool enableCBMC = false;           // Configurational bias MC
        int cbmcTrials = 10;               // Number of CBMC trials
        bool enableCavityBias = false;     // Cavity bias insertion
        double cavityRadius = 2.5;         // Cavity radius for bias
        
        // Adaptive sampling
        bool enableAdaptive = false;       // Adaptive move probabilities
        int adaptiveInterval = 1000;       // Update interval
        double adaptiveRate = 0.1;         // Learning rate
        
        // Multi-component options
        bool enableMultiComponent = false;  // Multi-component GCMC
        bool coupleInsertionDeletion = true; // Couple ins/del for charge neutrality
    };
    
    // Constructor and destructor
    explicit GCMCController(const Config& config = Config());
    ~GCMCController();
    
    // Initialize with components
    void initialize(model::montecarlo::MCState* state,
                   movement::gcmc::GCMCEngine* engine,
                   movement::MultiTypeReservoir* reservoir);
    
    // Move execution
    movement::MovementResult performMove();
    movement::MovementResult performInsertion();
    movement::MovementResult performDeletion();
    movement::MovementResult performTranslation();
    movement::MovementResult performRotation();
    
    // Move selection
    std::string selectMoveType();
    int selectFragmentType();
    int selectMolecule();
    
    // Move probabilities
    void setMoveProbabilities(double insertion, double deletion,
                             double translation, double rotation);
    std::map<std::string, double> getMoveProbabilities() const;
    void updateAdaptiveProbabilities();
    
    // Chemical potential management
    void setChemicalPotential(const std::string& fragmentType, double mu);
    void setChemicalPotentials(const std::map<std::string, double>& mus);
    std::map<std::string, double> getChemicalPotentials() const;
    
    // Fragment management
    void addFragmentType(const std::string& name,
                        const std::vector<model::montecarlo::MCAtom>& atoms,
                        double chemicalPotential);
    std::vector<std::string> getFragmentTypes() const;
    
    // Statistics
    const movement::gcmc::GCMCStatistics& getStatistics() const;
    void resetStatistics();
    
    // Acceptance criteria
    void setAcceptanceCriteria(movement::gcmc::GCMCAcceptance* acceptance);
    double calculateAcceptanceProbability(const movement::MovementResult& result) const;
    
    // Configuration
    const Config& getConfig() const { return config_; }
    void updateConfig(const Config& config);
    
    // Multi-component control
    void enableMultiComponent(bool enable);
    void setCoupledInsertion(bool coupled);
    std::pair<std::string, std::string> selectCoupledFragments();
    
private:
    Config config_;
    model::montecarlo::MCState* state_;
    movement::gcmc::GCMCEngine* engine_;
    movement::MultiTypeReservoir* reservoir_;
    movement::gcmc::GCMCAcceptance* acceptance_;
    
    bool initialized_;
    
    // Move probabilities
    std::map<std::string, double> moveProbabilities_;
    std::map<std::string, int> moveAttempts_;
    std::map<std::string, int> moveAccepted_;
    
    // Chemical potentials
    std::map<std::string, double> chemicalPotentials_;
    
    // Statistics
    std::unique_ptr<movement::gcmc::GCMCStatistics> statistics_;
    
    // Random number generation
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniformDist_;
    
    // Private methods
    void checkInitialized() const;
    void normalizeProbabilities();
    void updateMoveStatistics(const std::string& moveType, bool accepted);
    double computeAdaptiveProbability(const std::string& moveType);
};

} // namespace gcmc
} // namespace simulation
} // namespace cpu
} // namespace platform
} // namespace pygcmc