// modules/core/include/pygcmc/core/force_field.hpp

#ifndef PYGCMC_CORE_FORCE_FIELD_HPP
#define PYGCMC_CORE_FORCE_FIELD_HPP

#include <vector>
#include <array>
#include <memory>
#include <string>
#include <cmath>
#include <stdexcept>
#include "pygcmc/core/force.hpp"
#include "pygcmc/core/system.hpp"

namespace pygcmc {
namespace core {

/**
 * @brief Lennard-Jones force
 */
class LennardJonesForce : public Force {
public:
    LennardJonesForce(double epsilon = 1.0, double sigma = 1.0)
        : epsilon_(epsilon), sigma_(sigma) {
        if (epsilon < 0.0 || sigma < 0.0) {
            throw std::invalid_argument("Force field parameters must be non-negative");
        }
    }

    double calculate_forces(const System& system,
                          std::vector<std::vector<std::array<double, 3>>>& forces) const override;

    bool uses_periodic_boundary_conditions() const override {
        return true;
    }

private:
    double epsilon_;
    double sigma_;
};

/**
 * @brief Coulomb force
 */
class CoulombForce : public Force {
public:
    CoulombForce(double k = 1.0) : k_(k) {
        if (k < 0.0) {
            throw std::invalid_argument("Coulomb constant must be non-negative");
        }
    }

    double calculate_forces(const System& system,
                          std::vector<std::vector<std::array<double, 3>>>& forces) const override;

    bool uses_periodic_boundary_conditions() const override {
        return true;
    }

private:
    double k_;
};

/**
 * @brief Represents a bond between two residues
 */
struct ResidueBond {
    size_t residue1;
    size_t residue2;
    double k;
    double r0;

    ResidueBond(size_t r1, size_t r2, double k_, double r0_)
        : residue1(r1), residue2(r2), k(k_), r0(r0_) {}
};

/**
 * @brief Force between bonded residues
 */
class ResidueBondedForce : public Force {
public:
    ResidueBondedForce() = default;

    void add_bond(size_t residue1, size_t residue2, double k, double r0) {
        bonds_.emplace_back(residue1, residue2, k, r0);
    }

    double calculate_forces(const System& system,
                          std::vector<std::vector<std::array<double, 3>>>& forces) const override;

    bool uses_periodic_boundary_conditions() const override {
        return true;
    }

private:
    std::vector<ResidueBond> bonds_;
};

/**
 * @brief Represents an angle between three residues
 */
struct ResidueAngle {
    size_t residue1;
    size_t residue2;
    size_t residue3;
    double k;
    double theta0;

    ResidueAngle(size_t r1, size_t r2, size_t r3, double k_, double theta0_)
        : residue1(r1), residue2(r2), residue3(r3), k(k_), theta0(theta0_) {}
};

/**
 * @brief Force for angle between three residues
 */
class ResidueAngleForce : public Force {
public:
    ResidueAngleForce() = default;

    void add_angle(size_t residue1, size_t residue2, size_t residue3,
                  double k, double theta0) {
        angles_.emplace_back(residue1, residue2, residue3, k, theta0);
    }

    double calculate_forces(const System& system,
                          std::vector<std::vector<std::array<double, 3>>>& forces) const override;

    bool uses_periodic_boundary_conditions() const override {
        return true;
    }

private:
    std::vector<ResidueAngle> angles_;
};

} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_FORCE_FIELD_HPP