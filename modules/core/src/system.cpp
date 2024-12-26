// modules/core/src/system.cpp
#include "pygcmc/core/system.hpp"
#include <iostream> // 用于调试
#include <array>

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

            // New: Added Coulomb potential calculation
            double coulomb = (particles_[i].charge * particles_[j].charge) / distance;

            // Updated: Total energy now includes both potentials
            double total_pair_energy = lj + coulomb;
            total_energy += total_pair_energy;

            // Updated debug output to include Coulomb energy
            std::cout << "Pair (" << i << ", " << j << "): "
                      << "Distance = " << distance << ", "
                      << "inv_r6 = " << inv_r6 << ", "
                      << "inv_r12 = " << inv_r12 << ", "
                      << "LJ Energy = " << lj << ", "
                      << "Coulomb Energy = " << coulomb << ", "
                      << "Total Pair Energy = " << total_pair_energy << std::endl;
        }
    }
    std::cout << "Total Energy: " << total_energy << std::endl;
    return total_energy;
}

size_t System::get_particle_count() const {
    return particles_.size();
}

void System::update_positions(double dt) {
    for (auto& particle : particles_) {
        particle.x += particle.vx * dt;
        particle.y += particle.vy * dt;
        particle.z += particle.vz * dt;
    }
}

void System::update_velocities(double dt) {
    std::vector<std::array<double, 3>> accelerations(particles_.size(), {0.0, 0.0, 0.0});

    for (size_t i = 0; i < particles_.size(); ++i) {
        for (size_t j = i + 1; j < particles_.size(); ++j) {
            double dx = particles_[i].x - particles_[j].x;
            double dy = particles_[i].y - particles_[j].y;
            double dz = particles_[i].z - particles_[j].z;
            double distance = std::sqrt(dx * dx + dy * dy + dz * dz) + 1e-12;
            
            double inv_r = sigma_ / distance;
            double inv_r6 = std::pow(inv_r, 6);
            double inv_r12 = inv_r6 * inv_r6;

            // Calculate forces
            double force_lj = 24.0 * epsilon_ / distance * (2.0 * inv_r12 - inv_r6);
            double force_coulomb = (particles_[i].charge * particles_[j].charge) / (distance * distance);
            double total_force = force_lj + force_coulomb;

            // Update accelerations
            double fx = total_force * (dx / distance);
            double fy = total_force * (dy / distance);
            double fz = total_force * (dz / distance);

            accelerations[i][0] += fx;
            accelerations[i][1] += fy;
            accelerations[i][2] += fz;

            accelerations[j][0] -= fx;
            accelerations[j][1] -= fy;
            accelerations[j][2] -= fz;
        }
    }

    // Update velocities
    for (size_t i = 0; i < particles_.size(); ++i) {
        particles_[i].vx += accelerations[i][0] * dt;
        particles_[i].vy += accelerations[i][1] * dt;
        particles_[i].vz += accelerations[i][2] * dt;
    }
}
