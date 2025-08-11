#pragma once
#include "../DrudeInterface.hpp"
#include "../DrudeStructures.hpp"
#include <vector>
#include <array>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace exp {

// Algorithm selection for Drude model
enum class DrudeAlgorithm {
    S1_POINT_CHARGE,     // CHARMM/OpenMM standard: 4-point charge with S1 screening
    S3S5_DIPOLE_TENSOR,  // Dipole tensor with S3/S5 screening
    S1_DIPOLE_FIELD,     // Hybrid: dipole field with S1 screening
    DIRECT_COULOMB       // No screening (for testing)
};

class DrudeSCFOM {
public:
    // Set the algorithm to use
    void setAlgorithm(DrudeAlgorithm algo) { algorithm = algo; }
    DrudeAlgorithm getAlgorithm() const { return algorithm; }
    
    // Main optimization function
    bool optimize(model::MCState& state,
                  const std::vector<DrudeParticle>& particles,
                  const std::vector<ScreenedPair>& pairs,
                  const DrudeSCFParams& params);

    // Utility function for PBC
    static void applyPBC(double& dx, double& dy, double& dz, 
                         const std::array<double, 3>& box);

private:
    // Calculate total electric field at each Drude particle
    void calculateElectricField(const model::MCState& state,
                                const std::vector<DrudeParticle>& particles,
                                const std::vector<ScreenedPair>& pairs,
                                std::vector<Vec3>& electricField,
                                const DrudeSCFParams& params) const;

    // Calculate external field from non-Drude charges
    void calculateExternalField(const model::MCState& state,
                                const std::vector<DrudeParticle>& particles,
                                const std::vector<ScreenedPair>& pairs,
                                std::vector<Vec3>& electricField,
                                const DrudeSCFParams& params) const;

    // Calculate induced field using CHARMM/OpenMM point charge model with Thole screening
    void calculateInducedField(const model::MCState& state,
                               const std::vector<DrudeParticle>& particles,
                               const std::vector<ScreenedPair>& pairs,
                               std::vector<Vec3>& electricField) const;

    // Update Drude positions based on electric field
    double updateDrudePositions(model::MCState& state,
                                const std::vector<DrudeParticle>& particles,
                                const std::vector<Vec3>& electricField,
                                double damping,
                                double maxStep,
                                double hardWall) const;

    // Check if two atoms are in the same molecule
    bool inSameMolecule(int atom1, int atom2, const model::MCState& state) const;

    // Calculate Thole S1 screening function (CHARMM/OpenMM standard)
    double tholeS1(double r, double alpha_i, double alpha_j, double thole_sum) const;
    
    // Calculate Thole S3 screening function for 1/r^3 dipole-dipole interactions
    double tholeS3(double r, double alpha_i, double alpha_j, double thole_sum) const;
    
    // Calculate Thole S5 screening function for 1/r^5 tensor component
    double tholeS5(double r, double alpha_i, double alpha_j, double thole_sum) const;
    
    // Calculate induced field using S1 point charge model (CHARMM standard)
    void calculateInducedFieldS1(const model::MCState& state,
                                  const std::vector<DrudeParticle>& particles,
                                  const std::vector<ScreenedPair>& pairs,
                                  std::vector<Vec3>& electricField) const;
    
    // Calculate induced field using S3/S5 dipole tensor model
    void calculateInducedFieldS3S5(const model::MCState& state,
                                    const std::vector<DrudeParticle>& particles,
                                    const std::vector<ScreenedPair>& pairs,
                                    std::vector<Vec3>& electricField) const;

    // Calculate spring energy for convergence check
    double calculateSpringEnergy(const model::MCState& state,
                                 const std::vector<DrudeParticle>& particles) const;

    // Store best positions for fallback
    struct DrudePosition {
        double dx, dy, dz;
    };
    
    // DIIS acceleration data structures
    struct DIISData {
        std::vector<std::vector<double>> positions;  // History of positions
        std::vector<std::vector<double>> residuals;  // History of residuals
        int historySize = 0;
        int maxHistory = 5;
        bool enabled = false;
    };
    
    // DIIS mixing method
    bool applyDIIS(DIISData& diis, 
                   const std::vector<DrudeParticle>& particles,
                   const std::vector<Vec3>& electricField,
                   model::MCState& state) const;
                   
private:
    // Current algorithm selection
    DrudeAlgorithm algorithm = DrudeAlgorithm::S1_POINT_CHARGE;  // Default to CHARMM standard
};

} // namespace exp
} // namespace cpu
} // namespace platform
} // namespace pygcmc