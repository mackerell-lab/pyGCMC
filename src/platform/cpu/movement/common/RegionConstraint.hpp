#ifndef PYGCMC_PLATFORM_CPU_MOVEMENT_REGION_CONSTRAINT_HPP
#define PYGCMC_PLATFORM_CPU_MOVEMENT_REGION_CONSTRAINT_HPP

#include <string>
#include <vector>
#include <memory>
#include "MovementUtils.hpp"

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

/**
 * Base class for region constraints
 */
class RegionConstraint {
public:
    virtual ~RegionConstraint() = default;

    // Check if a position is within the constraint region
    virtual bool isInRegion(const Vector3& position) const = 0;

    // Sample a random position within the region
    virtual Vector3 samplePosition() const = 0;

    // Get the volume of the region (for bias calculations)
    virtual double getVolume() const = 0;

    // Parse region string and create appropriate constraint
    static std::unique_ptr<RegionConstraint> parseRegion(const std::string& regionSpec,
                                                         const Vector3& boxSize);
};

/**
 * Spherical region constraint
 */
class SphereRegion : public RegionConstraint {
public:
    SphereRegion(const Vector3& center, double radius, const Vector3& boxSize)
        : center_(center), radius_(radius), boxSize_(boxSize) {}

    bool isInRegion(const Vector3& position) const override;
    Vector3 samplePosition() const override;
    double getVolume() const override;

private:
    Vector3 center_;
    double radius_;
    Vector3 boxSize_;
};

/**
 * Box region constraint
 */
class BoxRegion : public RegionConstraint {
public:
    BoxRegion(const Vector3& min, const Vector3& max, const Vector3& boxSize)
        : min_(min), max_(max), boxSize_(boxSize) {}

    bool isInRegion(const Vector3& position) const override;
    Vector3 samplePosition() const override;
    double getVolume() const override;

private:
    Vector3 min_;
    Vector3 max_;
    Vector3 boxSize_;
};

/**
 * Cylindrical region constraint
 */
class CylinderRegion : public RegionConstraint {
public:
    CylinderRegion(const Vector3& center, double radius, double height,
                   int axis, const Vector3& boxSize)
        : center_(center), radius_(radius), height_(height),
          axis_(axis), boxSize_(boxSize) {}

    bool isInRegion(const Vector3& position) const override;
    Vector3 samplePosition() const override;
    double getVolume() const override;

private:
    Vector3 center_;
    double radius_;
    double height_;
    int axis_; // 0=x, 1=y, 2=z
    Vector3 boxSize_;
};

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc

#endif // PYGCMC_PLATFORM_CPU_MOVEMENT_REGION_CONSTRAINT_HPP