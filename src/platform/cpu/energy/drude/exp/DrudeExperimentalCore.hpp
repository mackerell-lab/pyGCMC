#pragma once
#include "../DrudeInterface.hpp"
#include "../DrudeStructures.hpp"
#include <vector>
#include <memory>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace exp {

class DrudeSCFOpenMM; // forward declaration

class DrudeExperimentalCore final : public DrudeInterface {
public:
    DrudeExperimentalCore();
    ~DrudeExperimentalCore() = default;

    // Main interface methods
    double calculateEnergy(model::MCState& state) override;
    void calculateForces(model::MCState& state, std::vector<Vec3>& forces) override;

    // Particle management
    int addParticle(const DrudeParticle& particle) override;
    void addScreenedPair(const ScreenedPair& pair) override;

    // Extra API for experimental path
    void autoScreenPairs(double thole, double cutoff_nm);
    
    // Configuration
    void setAlgorithm(DrudeAlgorithm algorithm) override; // kept for interface
    void setParameters(const DrudeSCFParams& params) override;
    void clear() override;
    size_t getNumParticles() const override { return m_particles.size(); }
    size_t getNumScreenedPairs() const { return m_screenedPairs.size(); }

    // Enable/disable Coulomb energy (for testing)
    void setIncludeCoulomb(bool include) { m_includeCoulomb = include; }
    bool getIncludeCoulomb() const { return m_includeCoulomb; }

private:
    friend class DrudeSCFOpenMM;

    std::vector<DrudeParticle> m_particles;
    std::vector<ScreenedPair>  m_screenedPairs;
    DrudeSCFParams m_params;
    bool m_includeCoulomb;

    // New SCF optimizer only
    std::unique_ptr<DrudeSCFOpenMM> m_scf;

    // Helper to check if atoms are in same molecule
    bool inSameMolecule(int atom1, int atom2, const model::MCState& state) const;
};

} // namespace exp
} // namespace cpu
} // namespace platform
} // namespace pygcmc