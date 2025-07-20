#pragma once

#include "model/ModelModule.hpp"
#include "platform/platform.hpp"
#include "../common/EnergyConstants.hpp"
#include "../common/EnergyUtils.hpp"
#include <vector>
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief Drude oscillator parameters for a single particle
 */
struct DrudeParticle {
    int drudeIndex;      // Index of the Drude particle
    int parentIndex;     // Index of the parent atom
    int aniso1Index;     // Index for first anisotropy axis (-1 if none)
    int aniso2Index;     // Index for second anisotropy axis (-1 if none)
    int aniso3Index;     // Index for third anisotropy axis (-1 if none)
    int aniso4Index;     // Index for fourth anisotropy axis (-1 if none)
    double charge;       // Charge on the Drude particle
    double polarizability; // Isotropic polarizability (α)
    double aniso12;      // Anisotropy scale factor for axis 1-2
    double aniso34;      // Anisotropy scale factor for axis 3-4
    double kIsotropic;   // Force constant for isotropic term
    double kAniso1;      // Force constant for first anisotropic term
    double kAniso2;      // Force constant for second anisotropic term
};

/**
 * @brief Screened pair interaction between dipoles
 */
struct ScreenedPair {
    int dipole1;         // Index of first dipole in DrudeParticle array
    int dipole2;         // Index of second dipole in DrudeParticle array
    double thole;        // Thole parameter for screening
};

/**
 * @brief SCF parameters for Drude optimization
 */
struct DrudeSCFParams {
    double tolerance = 1.0;      // Force tolerance for convergence (kJ/mol/nm) - same as OpenMM default
    int maxIterations = 50;      // Maximum SCF iterations
    double dampingFactor = 0.5;  // Damping for large forces
    double forceCutoff = 10.0;   // Force cutoff for damping (in units of tolerance)
    double maxDrudeDistance = 0.02; // Maximum Drude-parent distance (nm) - hard wall constraint
};

/**
 * @brief OPT3 coefficients for Drude model
 */
struct OPT3Coefficients {
    double c0 = 0.0;    // Zero-order weight (typically 0 for Drude)
    double c1 = 1.0;    // First-order weight (base response)
    double c2 = 0.0;    // Second-order weight (often 0 for simplicity)
    double c3 = 0.0;    // Third-order weight (often 0 for simplicity)
};

/**
 * @brief OPT4 coefficients for Drude model
 */
struct OPT4Coefficients {
    double c0 = 0.0;    // Zero-order weight (typically 0 for Drude)
    double c1 = 1.0;    // First-order weight (base response)
    double c2 = 0.0;    // Second-order weight
    double c3 = 0.0;    // Third-order weight
    double c4 = 0.0;    // Fourth-order weight
};

/**
 * @brief Algorithm selection for Drude optimization
 */
enum class DrudeAlgorithm {
    SCF,           // Standard Self-Consistent Field
    OPT3,          // 3rd order optimization
    OPT4,          // 4th order optimization
    AdaptiveOPT,   // Adaptive order selection
    HybridOPT,     // OPT prediction + SCF refinement
    SmartOPT3,     // OPT3 with smart r0 correction
    FBP,           // Force Balance Predictor
    ConjugateGradient  // Matrix-free Conjugate Gradient method
};

/**
 * @brief Calculate Drude oscillator energy and forces using SCF method
 * 
 * This implements the classical Drude oscillator model for polarizable force fields.
 * The SCF (Self-Consistent Field) method iteratively minimizes the positions of
 * Drude particles to find the ground state energy.
 */
class DrudeForce {
    friend class DrudeConjugateGradient;
public:
    /**
     * @brief Constructor
     */
    DrudeForce();
    
    /**
     * @brief Add a Drude particle
     * @return Index of the added particle
     */
    int addParticle(int drudeIndex, int parentIndex, 
                   int aniso1Index, int aniso2Index,
                   int aniso3Index, int aniso4Index,
                   double charge, double polarizability,
                   double aniso12, double aniso34);
    
    /**
     * @brief Add a screened pair interaction
     */
    void addScreenedPair(int dipole1, int dipole2, double thole);
    
    /**
     * @brief Set SCF parameters
     */
    void setSCFParameters(const DrudeSCFParams& params) { scfParams = params; }
    
    /**
     * @brief Simple 3D vector class for forces
     */
    struct Vec3 {
        double x, y, z;
        
        Vec3() : x(0), y(0), z(0) {}
        Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
        
        Vec3 operator+(const Vec3& v) const { return Vec3(x+v.x, y+v.y, z+v.z); }
        Vec3 operator-(const Vec3& v) const { return Vec3(x-v.x, y-v.y, z-v.z); }
        Vec3 operator*(double s) const { return Vec3(x*s, y*s, z*s); }
        Vec3& operator+=(const Vec3& v) { x+=v.x; y+=v.y; z+=v.z; return *this; }
        Vec3& operator-=(const Vec3& v) { x-=v.x; y-=v.y; z-=v.z; return *this; }
        
        double dot(const Vec3& v) const { return x*v.x + y*v.y + z*v.z; }
        double norm2() const { return x*x + y*y + z*z; }
        double norm() const { return std::sqrt(norm2()); }
        Vec3 normalized() const { double n = norm(); return n > 0 ? (*this)*(1.0/n) : Vec3(); }
    };
    
    /**
     * @brief Calculate energy using SCF method
     * @param state The molecular state
     * @return Total Drude energy
     */
    double calculateEnergySCF(model::MCState& state);
    
    /**
     * @brief Calculate energy using OPT3 method
     * @param state The molecular state
     * @return Total Drude energy
     */
    double calculateEnergyOPT3(model::MCState& state);
    
    /**
     * @brief Enable/disable OPT3 algorithm
     */
    void setUseOPT3(bool use) { useOPT3 = use; }
    bool getUseOPT3() const { return useOPT3; }
    
    /**
     * @brief Set algorithm
     */
    void setAlgorithm(DrudeAlgorithm algo) { algorithm = algo; }
    DrudeAlgorithm getAlgorithm() const { return algorithm; }
    
    /**
     * @brief Set OPT3 coefficients
     */
    void setOPT3Coefficients(double c0, double c1, double c2, double c3);
    
    /**
     * @brief Get current OPT3 coefficients
     */
    OPT3Coefficients getOPT3Coefficients() const { return opt3Coeffs; }
    
    /**
     * @brief Set OPT4 coefficients
     */
    void setOPT4Coefficients(double c0, double c1, double c2, double c3, double c4);
    
    /**
     * @brief Get current OPT4 coefficients
     */
    OPT4Coefficients getOPT4Coefficients() const { return opt4Coeffs; }
    
    /**
     * @brief Training data for OPT3 optimization
     */
    struct OPT3TrainingData {
        std::vector<Vec3> r0;        // Zero-order displacements
        std::vector<Vec3> r1;        // First-order displacements
        std::vector<Vec3> r2;        // Second-order displacements
        std::vector<Vec3> r3;        // Third-order displacements
        std::vector<Vec3> r_scf;     // SCF converged positions (ground truth)
        std::vector<Vec3> parentPos; // Parent positions for reference
    };
    
    /**
     * @brief Collect training data for OPT3 optimization
     * @param state The molecular state
     * @return Training data containing perturbation orders and SCF solution
     */
    OPT3TrainingData collectTrainingData(model::MCState& state);
    
    /**
     * @brief Enable/disable training mode
     */
    void enableTrainingMode(bool enable) { trainingMode = enable; }
    bool isTrainingMode() const { return trainingMode; }
    
    /**
     * @brief Calculate forces on all particles
     * @param state The molecular state
     * @param forces Output force array (must be pre-allocated)
     */
    void calculateForces(model::MCState& state, std::vector<Vec3>& forces);
    
    /**
     * @brief Get number of Drude particles
     */
    int getNumParticles() const { return particles.size(); }
    
    /**
     * @brief Get number of screened pairs
     */
    int getNumScreenedPairs() const { return screenedPairs.size(); }
    
private:
    std::vector<DrudeParticle> particles;
    std::vector<ScreenedPair> screenedPairs;
    DrudeSCFParams scfParams;
    OPT3Coefficients opt3Coeffs;  // OPT3 coefficients
    OPT4Coefficients opt4Coeffs;  // OPT4 coefficients
    DrudeAlgorithm algorithm = DrudeAlgorithm::SCF;  // Algorithm selection
    bool useOPT3 = false;  // Legacy: Use OPT3 algorithm instead of SCF
    bool trainingMode = false;  // Enable training data collection
    
    // OPT specific methods
    bool minimizeDrudePositionsWithOPT3(model::MCState& state);
    bool minimizeDrudePositionsWithOPT4(model::MCState& state);
    bool minimizeDrudePositionsAdaptive(model::MCState& state);
    bool minimizeDrudePositionsHybrid(model::MCState& state);
    bool minimizeDrudePositionsWithSmartOPT3(model::MCState& state);
    
    // Smart r0 methods
    void calculateSmartR0(model::MCState& state, std::vector<Vec3>& r0);
    void calculateStandardR0(model::MCState& state, std::vector<Vec3>& r0);
    double calculateLocalDensity(const model::MCState& state, const model::MCAtom& drude, int drudeIdx);
    bool isInHydrogenBondNetwork(const model::MCState& state, const DrudeParticle& particle);
    double calculateLocalChargeAsymmetry(const model::MCState& state, const model::MCAtom& drude, int drudeIdx);
    void updateDrudePositions(model::MCState& state, const std::vector<Vec3>& displacements);
    void calculateElectricFieldResponse(model::MCState& state, std::vector<Vec3>& response);
    Vec3 calculateElectricFieldAtPoint(const model::MCState& state, const model::MCAtom& point, int excludeDrudeIdx);
    double calculateDrudeForces(const model::MCState& state, std::vector<Vec3>& drudeForces);
    void performSingleSCFIteration(model::MCState& state, const std::vector<Vec3>& drudeForces);
    double calculateEnergyDirect(const model::MCState& state);
    void calculateElectricFieldAtDrudes(model::MCState& state, std::vector<Vec3>& electricField, bool includeDrudes);
    static void applyPBC(double& dx, double& dy, double& dz, double boxX, double boxY, double boxZ);
    
    // Force Balance Predictor methods
    bool minimizeDrudePositionsWithFBP(model::MCState& state);
    void calculateFixedElectricField(const model::MCState& state, std::vector<Vec3>& fixedField);
    void calculateDrudeForces(const model::MCState& state, std::vector<Vec3>& forces, const std::vector<Vec3>& fixedField);
    
    // Conjugate Gradient method
    bool minimizeDrudePositionsWithCG(model::MCState& state);
    double calculateEnergyFBP(model::MCState& state);
    
    /**
     * @brief Minimize Drude particle positions using SCF
     * @param state The molecular state
     * @param forces Current forces on all particles
     * @return true if converged, false otherwise
     */
    bool minimizeDrudePositions(model::MCState& state, std::vector<Vec3>& forces);
    
    /**
     * @brief Calculate harmonic restraint energy and forces
     */
    double calculateHarmonicEnergy(model::MCState& state, std::vector<Vec3>& forces);
    
    /**
     * @brief Calculate screened Coulomb interactions
     */
    double calculateScreenedCoulombEnergy(model::MCState& state, std::vector<Vec3>& forces);
    
    /**
     * @brief Calculate all Coulomb interactions (for SCF)
     */
    double calculateCoulombEnergy(model::MCState& state, std::vector<Vec3>& forces);
};

// Forward declaration of helper function used in OPT3/OPT4
void calculateElectricFieldAtDrudes(
    const model::MCState& state,
    const std::vector<DrudeParticle>& particles,
    const std::vector<ScreenedPair>& screenedPairs,
    const std::vector<DrudeForce::Vec3>& drudePositions,
    std::vector<DrudeForce::Vec3>& electricField,
    bool includeStaticOnly = false);

} // namespace cpu
} // namespace platform
} // namespace pygcmc