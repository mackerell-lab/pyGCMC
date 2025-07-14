#ifndef PGPCONTEXT_HPP
#define PGPCONTEXT_HPP

#include <vector>
#include <complex>
#include <memory>
#include <cmath>
#include <array>
#include "model/montecarlo/MCMain.hpp"

// Forward declarations
namespace pygcmc {
namespace platform {
namespace cpu {
    struct PMEParams;
}
}
}

namespace platform {
namespace cpu {
namespace energy {
namespace pgp {

/**
 * PGPContext - Instance-based PGP implementation without any global state
 * 
 * This class encapsulates all PGP parameters and state, allowing multiple
 * independent instances to coexist without interference. This design eliminates
 * all race conditions and memory corruption issues in parallel execution.
 */
class PGPContext {
public:
    // Constructor
    PGPContext() = default;
    
    // Destructor
    ~PGPContext() = default;
    
    // Delete copy constructor and assignment to prevent accidental sharing
    PGPContext(const PGPContext&) = delete;
    PGPContext& operator=(const PGPContext&) = delete;
    
    // Allow move semantics
    PGPContext(PGPContext&&) = default;
    PGPContext& operator=(PGPContext&&) = default;
    
    /**
     * Initialize PGP parameters for this instance
     */
    void initialize(double cutoff, 
                   const std::array<double, 3>& box,
                   double alpha,
                   const std::array<int, 3>& meshSize,
                   double potential_cutoff,
                   const std::array<int, 3>& potentialGridSize,
                   int splineOrder,
                   double tolerance);
    
    /**
     * Compute system energy using this PGP instance
     */
    struct EnergyComponents {
        double total;
        double real_space;
        double reciprocal;
        double self;
        double vdw;
    };
    
    EnergyComponents computeSystemEnergy(pygcmc::model::montecarlo::MCState& state) const;
    
    /**
     * Compute movement energy using this PGP instance
     */
    EnergyComponents computeMovementEnergy(pygcmc::model::montecarlo::MCState& state, 
                                         const std::vector<int>& movementResidues) const;
    
    /**
     * Precompute grid potential (placeholder for future implementation)
     */
    void precomputeGridPotential(pygcmc::model::montecarlo::MCState& state, int atomType);
    
private:
    // PGP-specific parameters (no inheritance from PME!)
    struct Parameters {
        // Basic parameters
        double alpha{1.0};
        double cutoff{1.0};
        double potential_cutoff{1.0};
        std::array<double, 3> box{10.0, 10.0, 10.0};
        std::array<int, 3> meshSize{64, 64, 64};
        std::array<int, 3> potentialGridSize{64, 64, 64};
        int splineOrder{4};
        double tolerance{1e-6};
        
        // Derived parameters
        double alphaEwald{0.0};
        double kmax{0.0};
        double ewaldCoeff{0.0};
        double selfEnergyCoeff{0.0};
        double epsilon_r{1.0};
        
        // Instance-specific lookup tables (not shared!)
        std::vector<double> erfcTable;
        std::vector<double> expTable;
        std::vector<double> dExpTable;
        
        // Instance-specific grid structures
        std::vector<std::complex<double>> pmeGrid;
        std::vector<std::complex<double>> pmeGridSaved;
        std::vector<double> potentialGrid;
        std::vector<double> bsplineCoeffs;
        
        // Grid dimensions
        int gridSizeX{0};
        int gridSizeY{0};
        int gridSizeZ{0};
        int gridSizeYZ{0};
        int gridTotal{0};
        
        // Table dimensions
        static constexpr int ERFC_TABLE_SIZE = 2000;
        static constexpr int EXP_TABLE_SIZE = 2000;
        static constexpr double ERFC_TABLE_SCALE = 2.0;
        static constexpr double EXP_TABLE_SCALE = 10.0;
    } params_;
    
    // Private helper methods
    void setupTables();
    void setupGrids();
    double computeRealSpaceEnergy(pygcmc::model::montecarlo::MCState& state, bool movement_only) const;
    double computeReciprocalEnergy(pygcmc::model::montecarlo::MCState& state, bool movement_only) const;
    double computeSelfEnergy(pygcmc::model::montecarlo::MCState& state, bool movement_only) const;
    double computeVdWEnergy(pygcmc::model::montecarlo::MCState& state, bool movement_only) const;
    
    // Helper functions for table lookups
    double getErfcValue(double x) const;
    double getExpValue(double x) const;
    
    // New helper functions for FFT-based precomputation
    void precomputeGridPotentialSimple(pygcmc::model::montecarlo::MCState& state, int atomType);
    static void applyReciprocalSpacePotentialFactors(::pygcmc::platform::cpu::PMEParams& pme, const Parameters& params);
};

} // namespace pgp
} // namespace energy
} // namespace cpu
} // namespace platform

#endif // PGPCONTEXT_HPP