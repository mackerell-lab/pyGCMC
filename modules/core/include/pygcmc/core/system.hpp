// modules/core/include/pygcmc/core/system.hpp
#ifndef PYGCMC_CORE_SYSTEM_HPP
#define PYGCMC_CORE_SYSTEM_HPP

#include <vector>
#include <cmath>

// 粒子结构体
struct Particle {
    double x, y, z;    // 位置坐标
    int type;          // 粒子类型
    double charge;     // 粒子电荷

    Particle(double x_val = 0.0, double y_val = 0.0, double z_val = 0.0, int t = 0, double q = 0.0)
        : x(x_val), y(y_val), z(z_val), type(t), charge(q) {}
};

// 系统类管理整个粒子系统
class System {
public:
    System(double epsilon = 1.0, double sigma = 1.0); // 添加参数
    ~System();

    // 添加粒子
    void add_particle(const Particle& particle);

    // 移除粒子（按索引）
    void remove_particle(int index);

    // 计算总能量
    double compute_total_energy() const;

    // 获取粒子数量
    size_t get_particle_count() const;

private:
    std::vector<Particle> particles_;
    double epsilon_; // 势能参数
    double sigma_;   // 势能参数
};

#endif // PYGCMC_CORE_SYSTEM_HPP
