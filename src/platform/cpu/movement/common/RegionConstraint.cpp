#include "RegionConstraint.hpp"
#include <sstream>
#include <cmath>
#include <stdexcept>

namespace pygcmc {
namespace platform {
namespace cpu {
namespace movement {

// Parse region specification string
// Formats:
// - "sphere x y z r" - sphere centered at (x,y,z) with radius r (in nm)
// - "box x1 y1 z1 x2 y2 z2" - box from (x1,y1,z1) to (x2,y2,z2) (in nm)
// - "cylinder x y z r h axis" - cylinder centered at (x,y,z) with radius r, height h, along axis (x/y/z) (in nm)
// - "" or "all" - entire box (default)
std::unique_ptr<RegionConstraint> RegionConstraint::parseRegion(
    const std::string& regionSpec, const Vector3& boxSize) {

    if (regionSpec.empty() || regionSpec == "all") {
        // No constraint - entire box
        return std::make_unique<BoxRegion>(
            Vector3(0, 0, 0),
            boxSize,
            boxSize
        );
    }

    std::istringstream iss(regionSpec);
    std::string type;
    iss >> type;

    if (type == "sphere") {
        double x, y, z, r;
        if (!(iss >> x >> y >> z >> r)) {
            throw std::runtime_error("Invalid sphere region format. Expected: sphere x y z r");
        }
        // Input is already in nm (consistent with box_size, cutoff, etc.)
        Vector3 center(x, y, z);
        double radius = r;
        return std::make_unique<SphereRegion>(center, radius, boxSize);

    } else if (type == "box") {
        double x1, y1, z1, x2, y2, z2;
        if (!(iss >> x1 >> y1 >> z1 >> x2 >> y2 >> z2)) {
            throw std::runtime_error("Invalid box region format. Expected: box x1 y1 z1 x2 y2 z2");
        }
        // Input is already in nm (consistent with box_size, cutoff, etc.)
        Vector3 min(x1, y1, z1);
        Vector3 max(x2, y2, z2);
        return std::make_unique<BoxRegion>(min, max, boxSize);

    } else if (type == "cylinder") {
        double x, y, z, r, h;
        std::string axisStr;
        if (!(iss >> x >> y >> z >> r >> h >> axisStr)) {
            throw std::runtime_error("Invalid cylinder region format. Expected: cylinder x y z r h axis");
        }
        // Input is already in nm (consistent with box_size, cutoff, etc.)
        Vector3 center(x, y, z);
        double radius = r;
        double height = h;

        int axis = 2; // default z-axis
        if (axisStr == "x") axis = 0;
        else if (axisStr == "y") axis = 1;
        else if (axisStr == "z") axis = 2;
        else throw std::runtime_error("Invalid axis. Must be x, y, or z");

        return std::make_unique<CylinderRegion>(center, radius, height, axis, boxSize);

    } else {
        throw std::runtime_error("Unknown region type: " + type + ". Supported: sphere, box, cylinder");
    }
}

// SphereRegion implementation
bool SphereRegion::isInRegion(const Vector3& position) const {
    // Calculate distance from center with PBC
    Vector3 diff = position - center_;

    // Apply minimum image convention
    if (diff.x > boxSize_.x * 0.5) diff.x -= boxSize_.x;
    else if (diff.x < -boxSize_.x * 0.5) diff.x += boxSize_.x;

    if (diff.y > boxSize_.y * 0.5) diff.y -= boxSize_.y;
    else if (diff.y < -boxSize_.y * 0.5) diff.y += boxSize_.y;

    if (diff.z > boxSize_.z * 0.5) diff.z -= boxSize_.z;
    else if (diff.z < -boxSize_.z * 0.5) diff.z += boxSize_.z;

    double distSq = diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
    return distSq <= radius_ * radius_;
}

Vector3 SphereRegion::samplePosition() const {
    // Sample uniformly within sphere using rejection sampling
    while (true) {
        // Generate random position in cube
        Vector3 offset(
            utils::RandomUtils::uniform(-radius_, radius_),
            utils::RandomUtils::uniform(-radius_, radius_),
            utils::RandomUtils::uniform(-radius_, radius_)
        );

        // Check if within sphere
        double distSq = offset.x * offset.x + offset.y * offset.y + offset.z * offset.z;
        if (distSq <= radius_ * radius_) {
            Vector3 pos = center_ + offset;
            // Apply PBC to ensure position is within box
            while (pos.x < 0) pos.x += boxSize_.x;
            while (pos.x >= boxSize_.x) pos.x -= boxSize_.x;
            while (pos.y < 0) pos.y += boxSize_.y;
            while (pos.y >= boxSize_.y) pos.y -= boxSize_.y;
            while (pos.z < 0) pos.z += boxSize_.z;
            while (pos.z >= boxSize_.z) pos.z -= boxSize_.z;
            return pos;
        }
    }
}

double SphereRegion::getVolume() const {
    return (4.0 / 3.0) * M_PI * radius_ * radius_ * radius_;
}

// BoxRegion implementation
bool BoxRegion::isInRegion(const Vector3& position) const {
    // Check if position is within box bounds
    return position.x >= min_.x && position.x <= max_.x &&
           position.y >= min_.y && position.y <= max_.y &&
           position.z >= min_.z && position.z <= max_.z;
}

Vector3 BoxRegion::samplePosition() const {
    return Vector3(
        utils::RandomUtils::uniform(min_.x, max_.x),
        utils::RandomUtils::uniform(min_.y, max_.y),
        utils::RandomUtils::uniform(min_.z, max_.z)
    );
}

double BoxRegion::getVolume() const {
    return (max_.x - min_.x) * (max_.y - min_.y) * (max_.z - min_.z);
}

// CylinderRegion implementation
bool CylinderRegion::isInRegion(const Vector3& position) const {
    // Check height constraint along axis
    double axialPos, axialMin, axialMax;
    if (axis_ == 0) {
        axialPos = position.x;
        axialMin = center_.x - height_ * 0.5;
        axialMax = center_.x + height_ * 0.5;
    } else if (axis_ == 1) {
        axialPos = position.y;
        axialMin = center_.y - height_ * 0.5;
        axialMax = center_.y + height_ * 0.5;
    } else {
        axialPos = position.z;
        axialMin = center_.z - height_ * 0.5;
        axialMax = center_.z + height_ * 0.5;
    }

    if (axialPos < axialMin || axialPos > axialMax) {
        return false;
    }

    // Check radial constraint in perpendicular plane
    Vector3 diff = position - center_;

    // Apply minimum image convention for radial components
    if (axis_ != 0) {
        if (diff.x > boxSize_.x * 0.5) diff.x -= boxSize_.x;
        else if (diff.x < -boxSize_.x * 0.5) diff.x += boxSize_.x;
    }
    if (axis_ != 1) {
        if (diff.y > boxSize_.y * 0.5) diff.y -= boxSize_.y;
        else if (diff.y < -boxSize_.y * 0.5) diff.y += boxSize_.y;
    }
    if (axis_ != 2) {
        if (diff.z > boxSize_.z * 0.5) diff.z -= boxSize_.z;
        else if (diff.z < -boxSize_.z * 0.5) diff.z += boxSize_.z;
    }

    double radialDistSq = 0.0;
    if (axis_ != 0) radialDistSq += diff.x * diff.x;
    if (axis_ != 1) radialDistSq += diff.y * diff.y;
    if (axis_ != 2) radialDistSq += diff.z * diff.z;

    return radialDistSq <= radius_ * radius_;
}

Vector3 CylinderRegion::samplePosition() const {
    // Sample uniformly within cylinder

    // Sample along axis
    double axialPos;
    if (axis_ == 0) {
        axialPos = utils::RandomUtils::uniform(
            center_.x - height_ * 0.5,
            center_.x + height_ * 0.5
        );
    } else if (axis_ == 1) {
        axialPos = utils::RandomUtils::uniform(
            center_.y - height_ * 0.5,
            center_.y + height_ * 0.5
        );
    } else {
        axialPos = utils::RandomUtils::uniform(
            center_.z - height_ * 0.5,
            center_.z + height_ * 0.5
        );
    }

    // Sample in radial plane using rejection sampling
    Vector3 pos = center_;
    if (axis_ == 0) pos.x = axialPos;
    else if (axis_ == 1) pos.y = axialPos;
    else pos.z = axialPos;

    while (true) {
        // Generate random position in square
        Vector3 radialOffset(0, 0, 0);
        if (axis_ != 0) radialOffset.x = utils::RandomUtils::uniform(-radius_, radius_);
        if (axis_ != 1) radialOffset.y = utils::RandomUtils::uniform(-radius_, radius_);
        if (axis_ != 2) radialOffset.z = utils::RandomUtils::uniform(-radius_, radius_);

        // Check if within circle
        double radialDistSq = 0.0;
        if (axis_ != 0) radialDistSq += radialOffset.x * radialOffset.x;
        if (axis_ != 1) radialDistSq += radialOffset.y * radialOffset.y;
        if (axis_ != 2) radialDistSq += radialOffset.z * radialOffset.z;

        if (radialDistSq <= radius_ * radius_) {
            pos = pos + radialOffset;
            // Apply PBC to ensure position is within box
            while (pos.x < 0) pos.x += boxSize_.x;
            while (pos.x >= boxSize_.x) pos.x -= boxSize_.x;
            while (pos.y < 0) pos.y += boxSize_.y;
            while (pos.y >= boxSize_.y) pos.y -= boxSize_.y;
            while (pos.z < 0) pos.z += boxSize_.z;
            while (pos.z >= boxSize_.z) pos.z -= boxSize_.z;
            return pos;
        }
    }
}

double CylinderRegion::getVolume() const {
    return M_PI * radius_ * radius_ * height_;
}

} // namespace movement
} // namespace cpu
} // namespace platform
} // namespace pygcmc