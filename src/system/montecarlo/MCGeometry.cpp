#include "MCGeometry.hpp"
#include <cmath>

namespace pygcmc {
namespace system {
namespace montecarlo {

float MCGeometry::getMinImageDistSqr(const model::MCState& state, float dx, float dy, float dz) const {
    // All distances are in nm, no conversion needed
    dx -= state.info.box[0] * std::round(dx / state.info.box[0]);
    dy -= state.info.box[1] * std::round(dy / state.info.box[1]);
    dz -= state.info.box[2] * std::round(dz / state.info.box[2]);
    return dx*dx + dy*dy + dz*dz;  // Returns square of distance in nm²
}

void MCGeometry::applyPBC(const model::MCState& state, float& x, float& y, float& z) const {
    // Box dimensions and coordinates are in nm, no conversion needed
    x -= state.info.box[0] * std::floor(x / state.info.box[0]);
    y -= state.info.box[1] * std::floor(y / state.info.box[1]);
    z -= state.info.box[2] * std::floor(z / state.info.box[2]);
}

void MCGeometry::updateGeometricCenter(const model::MCState& state, model::MCResidue& res) const {
    res.center[0] = res.center[1] = res.center[2] = 0.0f;

    // All coordinates are already in nm, no conversion needed
    for (int i = 0; i < res.atomCount; ++i) {
        const model::MCAtom& atom = state.atoms[res.atomStart + i];
        res.center[0] += atom.x;
        res.center[1] += atom.y;
        res.center[2] += atom.z;
    }

    if (res.atomCount > 0) {
        float invCount = 1.0f / res.atomCount;
        res.center[0] *= invCount;
        res.center[1] *= invCount;
        res.center[2] *= invCount;
    }
}

float MCGeometry::getDistance(const model::MCState& state,
                             float x1, float y1, float z1,
                             float x2, float y2, float z2) const {
    float dx = x2 - x1;
    float dy = y2 - y1;
    float dz = z2 - z1;

    return std::sqrt(getMinImageDistSqr(state, dx, dy, dz));
}

float MCGeometry::getDistanceSquared(const model::MCState& state,
                                    float x1, float y1, float z1,
                                    float x2, float y2, float z2) const {
    float dx = x2 - x1;
    float dy = y2 - y1;
    float dz = z2 - z1;

    return getMinImageDistSqr(state, dx, dy, dz);
}

} // namespace montecarlo
} // namespace system
} // namespace pygcmc
