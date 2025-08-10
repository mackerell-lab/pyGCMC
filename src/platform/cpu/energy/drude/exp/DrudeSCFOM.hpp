#pragma once
#include "../DrudeInterface.hpp"
#include "../DrudeStructures.hpp"
#include <vector>
#include <array>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace exp {

class DrudeSCFOM {
public:
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

    // Calculate Thole S3 screening function for 1/r^3 dipole-dipole interactions
    double tholeS3(double r, double alpha_i, double alpha_j, double thole_sum) const;
    
    // Calculate Thole S5 screening function for 1/r^5 tensor component
    double tholeS5(double r, double alpha_i, double alpha_j, double thole_sum) const;

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
};

} // namespace exp
} // namespace cpu
} // namespace platform
} // namespace pygcmc