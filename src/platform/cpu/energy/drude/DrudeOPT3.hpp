#pragma once

/**
 * @file DrudeOPT3.hpp
 * @brief Third-order perturbation theory optimizer for Drude oscillators
 */

#include "DrudeInterface.hpp"
#include "DrudeStructures.hpp"
#include <vector>

namespace pygcmc {
namespace platform {
namespace cpu {

/**
 * @brief OPT3 optimizer for Drude positions
 *
 * Uses perturbation theory expansion to approximate Drude positions:
 * r = c0*r0 + c1*r1 + c2*r2 + c3*r3
 */
class DrudeOPT3 : public DrudeOptimizer {
public:
    DrudeOPT3() = default;
    ~DrudeOPT3() = default;

    bool optimize(
        model::MCState& state,
        const std::vector<DrudeParticle>& particles,
        const std::vector<ScreenedPair>& screenedPairs,
        const DrudeSCFParams& params
    ) override;

    const char* getName() const override { return "OPT3"; }

    /**
     * @brief Set expansion coefficients
     */
    void setCoefficients(const OPT3Coefficients& coeffs) {
        m_coefficients = coeffs;
    }

private:
    OPT3Coefficients m_coefficients;

    /**
     * @brief Calculate electric field at atom positions
     * @param state MCState containing atom positions and charges
     * @param fields Output array for electric fields
     * @param particles Drude particles
     * @param screenedPairs Pairs with Thole screening
     * @param includeDrudes Whether to include fields from Drude particles
     */
    void calculateElectricField(
        const model::MCState& state,
        std::vector<Vec3>& fields,
        const std::vector<DrudeParticle>& particles,
        const std::vector<ScreenedPair>& screenedPairs,
        bool includeDrudes
    ) const;

    /**
     * @brief Check if two atoms are in the same molecule/residue
     */
    bool inSameMolecule(int atom1, int atom2, const model::MCState& state) const;
};

} // namespace cpu
} // namespace platform
} // namespace pygcmc
