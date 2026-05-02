/**
 * @file DrudeDirectPolarization.hpp
 * @brief Direct polarization approximation for Drude oscillators
 *
 * Ignores induced-induced interactions for maximum speed.
 * Suitable for GCMC screening where 5-10% error is acceptable.
 */

#ifndef PYGCMC_DRUDE_DIRECT_POLARIZATION_HPP
#define PYGCMC_DRUDE_DIRECT_POLARIZATION_HPP

#include "DrudeInterface.hpp"
#include "DrudeStructures.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {

class DrudeDirectPolarization : public DrudeOptimizer {
public:
    DrudeDirectPolarization() = default;
    ~DrudeDirectPolarization() override = default;

    bool optimize(
        model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<ScreenedPair>& screenedPairs,
        const DrudeSCFParams& params
    ) override;

    const char* getName() const override { return "DirectPolarization"; }

private:
    /**
     * Compute permanent electric field at Drude parent positions
     * Only includes fixed charges (non-Drude atoms)
     */
    void computePermanentFields(
        const model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        std::vector<Vec3>& fields
    );

    /**
     * Check if two atoms are in the same molecule
     */
    bool inSameMolecule(int atom1, int atom2, const model::MCState& state);
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_DRUDE_DIRECT_POLARIZATION_HPP
