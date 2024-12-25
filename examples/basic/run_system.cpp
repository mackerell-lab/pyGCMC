// examples/basic/run_system.cpp
#include <iostream>
#include "pygcmc/core/system.hpp"

int main() {
    // 使用默认的 epsilon=1.0 和 sigma=1.0
    System system;

    // 添加粒子，使得距离约为 3.355
    system.add_particle(Particle(0.0, 0.0, 0.0, 1, 0.5));
    system.add_particle(Particle(3.355, 0.0, 0.0, 2, -0.5));

    // 输出粒子数量
    std::cout << "Particle count: " << system.get_particle_count() << std::endl;

    // 计算并输出总能量
    double energy = system.compute_total_energy();
    std::cout << "Total energy: " << energy << std::endl;

    return 0;
}
