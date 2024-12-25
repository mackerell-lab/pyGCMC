// modules/core/src/system.cpp
#include "pygcmc/core/system.hpp"
#include <iostream> // 用于调试

System::System(double epsilon, double sigma)
    : epsilon_(epsilon), sigma_(sigma) {
    // 初始化系统，如果需要的话
}

System::~System() {
    // 清理资源，如果需要的话
}

void System::add_particle(const Particle& particle) {
    particles_.emplace_back(particle);
}

void System::remove_particle(int index) {
    if (index >= 0 && index < static_cast<int>(particles_.size())) {
        particles_.erase(particles_.begin() + index);
    }
}

double System::compute_total_energy() const {
    double total_energy = 0.0;
    // Lennard-Jones 势能计算
    for (size_t i = 0; i < particles_.size(); ++i) {
        for (size_t j = i + 1; j < particles_.size(); ++j) {
            double dx = particles_[i].x - particles_[j].x;
            double dy = particles_[i].y - particles_[j].y;
            double dz = particles_[i].z - particles_[j].z;
            double distance = std::sqrt(dx * dx + dy * dy + dz * dz) + 1e-12; // 避免除零

            // 计算 (sigma / r)^6 和 (sigma / r)^12
            double inv_r = sigma_ / distance;
            double inv_r6 = std::pow(inv_r, 6);
            double inv_r12 = inv_r6 * inv_r6;

            // Lennard-Jones 势能
            double lj = 4.0 * epsilon_ * (inv_r12 - inv_r6);
            total_energy += lj;

            // 调试输出
            std::cout << "Pair (" << i << ", " << j << "): "
                      << "Distance = " << distance << ", "
                      << "inv_r6 = " << inv_r6 << ", "
                      << "inv_r12 = " << inv_r12 << ", "
                      << "LJ Energy = " << lj << std::endl;
        }
    }
    std::cout << "Total Energy: " << total_energy << std::endl;
    return total_energy;
}

size_t System::get_particle_count() const {
    return particles_.size();
}
