#include "energyLJ.hpp"
#include <cmath>

namespace pygcmc {
namespace platform {
namespace cpu {

// 从energyCommon.hpp中移植过来的实现
float calculateSwitchingFunction(float r, const model::MCInfo& info) {
    if (!info.use_switching || r <= info.r_on) {
        return 1.0f;  // No switching below r_on
    }
    if (r >= info.r_off) {
        return 0.0f;  // Zero potential beyond r_off
    }
    
    // Calculate CHARMM-style switching function
    // S(r) = [(r_off^2 - r^2)^2 * (r_off^2 + 2r^2 - 3r_on^2)] / (r_off^2 - r_on^2)^3
    float r2 = r * r;
    float ron2 = info.r_on * info.r_on;
    float roff2 = info.r_off * info.r_off;
    
    float numerator = (roff2 - r2) * (roff2 - r2) * (roff2 + 2.0f*r2 - 3.0f*ron2);
    float denominator = (roff2 - ron2) * (roff2 - ron2) * (roff2 - ron2);
    
    return numerator / denominator;
}

// 特化的函数重载，用于处理混合类型参数的情况 - 实现
double calculateLJEnergy(
    double r2, 
    double sigma, 
    double eps, 
    const model::MCInfo& info,
    float min_safe_distance,
    float max_safe_energy
) {
    return calculateLJEnergy(r2, sigma, eps, info, 
                           static_cast<double>(min_safe_distance), 
                           static_cast<double>(max_safe_energy));
}

// 特化的函数重载，用于处理混合类型参数的情况 - 实现
float calculateLJEnergy(
    float r2, 
    float sigma, 
    float eps, 
    const model::MCInfo& info,
    double min_safe_distance,
    double max_safe_energy
) {
    return calculateLJEnergy(r2, sigma, eps, info, 
                           static_cast<float>(min_safe_distance), 
                           static_cast<float>(max_safe_energy));
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc 