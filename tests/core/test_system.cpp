// tests/core/test_system.cpp
#include <gtest/gtest.h>
#include "pygcmc/core/system.hpp"

namespace {

using namespace pygcmc::core;

class SystemTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up default system with unit parameters
        system = std::make_unique<System>(1.0, 1.0);
    }

    // Helper method to create a test particle
    Particle create_test_particle(int serial, const std::string& name, 
                                double x, double y, double z, 
                                double charge = 0.0) {
        // 确保所有必需字段都有有效值
        Particle p;
        p.serial = serial;
        p.name = name;
        p.residue = "TEST";  // 确保有残基名
        p.sequence = 1;      // 确保有序列号
        p.x = x;
        p.y = y;
        p.z = z;
        p.charge = charge;
        p.type = "TEST";     // 确保有类型
        return p;
    }

    std::unique_ptr<System> system;

    // Add new helper method for position/velocity verification
    void verify_position_velocity(const Particle& p, 
                                const std::array<double, 3>& expected_pos,
                                const std::array<double, 3>& expected_vel) {
        auto pos = p.position();
        auto vel = p.velocity();
        
        for (int i = 0; i < 3; ++i) {
            EXPECT_DOUBLE_EQ(pos[i], expected_pos[i]) 
                << "Position mismatch at index " << i;
            EXPECT_DOUBLE_EQ(vel[i], expected_vel[i]) 
                << "Velocity mismatch at index " << i;
        }
    }

    // Add helper method for center of mass verification
    void verify_center_of_mass(const Residue& res, 
                             const std::array<double, 3>& expected_com) {
        auto com = res.center_of_mass();
        for (int i = 0; i < 3; ++i) {
            EXPECT_DOUBLE_EQ(com[i], expected_com[i]) 
                << "Center of mass mismatch at index " << i;
        }
    }

    // Helper method to create a test residue
    Residue create_test_residue(const std::vector<Particle>& particles, 
                               const std::string& name, int seq_num, char chain) {
        // 确保序列号是正数
        if (seq_num <= 0) seq_num = 1;
        
        // 确保名称不为空
        std::string res_name = name.empty() ? "TEST" : name;
        
        Residue res(res_name, seq_num, chain);
        
        // 确保粒子的残基信息与残基匹配
        std::vector<Particle> valid_particles;
        for (auto p : particles) {
            p.residue = res_name;
            p.sequence = seq_num;
            if (p.is_valid()) {
                valid_particles.push_back(p);
            }
        }
        
        if (valid_particles.empty()) {
            throw std::runtime_error("No valid particles provided for residue");
        }
        
        res.atoms = valid_particles;
        return res;
    }
};

TEST_F(SystemTest, InitialState) {
    EXPECT_EQ(system->get_particle_count(), 0);
    
    auto [kinetic, potential] = system->get_system_state();
    EXPECT_DOUBLE_EQ(kinetic, 0.0);
    EXPECT_DOUBLE_EQ(potential, 0.0);
}

TEST_F(SystemTest, ParticleManagement) {
    // Add particle
    auto p1 = create_test_particle(1, "H1", 0.0, 0.0, 0.0, 0.5);
    system->add_particle(p1);
    EXPECT_EQ(system->get_particle_count(), 1);
    
    // Add another particle
    auto p2 = create_test_particle(2, "O1", 1.0, 0.0, 0.0, -1.0);
    system->add_particle(p2);
    EXPECT_EQ(system->get_particle_count(), 2);
    
    // Remove particle
    system->remove_particle(0);
    EXPECT_EQ(system->get_particle_count(), 1);
    
    // Try to remove invalid index
    EXPECT_NO_THROW(system->remove_particle(10));
    EXPECT_EQ(system->get_particle_count(), 1);
}

TEST_F(SystemTest, EnergyComputation) {
    // Add two particles at unit distance
    auto p1 = create_test_particle(1, "H1", 0.0, 0.0, 0.0, 0.5);
    auto p2 = create_test_particle(2, "O1", 1.0, 0.0, 0.0, -1.0);
    
    system->add_particle(p1);
    system->add_particle(p2);
    
    double energy = system->compute_total_energy();
    EXPECT_TRUE(std::isfinite(energy));
    EXPECT_NE(energy, 0.0);
}

TEST_F(SystemTest, PeriodicBoundary) {
    // Set up periodic boundary
    double box_size = 10.0;
    system->set_periodic_boundary(box_size);
    
    // Test PBC application
    EXPECT_DOUBLE_EQ(system->apply_pbc(11.0), 1.0);
    EXPECT_DOUBLE_EQ(system->apply_pbc(-1.0), 9.0);
    EXPECT_DOUBLE_EQ(system->apply_pbc(5.0), 5.0);
}

TEST_F(SystemTest, ParticleAccess) {
    auto p1 = create_test_particle(1, "H1", 0.0, 0.0, 0.0);
    system->add_particle(p1);
    
    // Test const access
    EXPECT_NO_THROW({
        const auto& particle = system->get_particle(0);
        EXPECT_EQ(particle.serial, 1);
        EXPECT_EQ(particle.name, "H1");
    });
    
    // Test non-const access
    EXPECT_NO_THROW({
        auto& particle = system->get_particle(0);
        particle.charge = 1.0;
        EXPECT_DOUBLE_EQ(particle.charge, 1.0);
    });
    
    // Test invalid access
    EXPECT_THROW(system->get_particle(1), std::out_of_range);
}

TEST_F(SystemTest, ParticlePositionVelocity) {
    // Create a particle with initial position and velocity
    Particle p1 = create_test_particle(1, "H1", 1.0, 2.0, 3.0);
    p1.set_velocity({0.5, -0.5, 1.0});

    // Test position method
    auto pos = p1.position();
    EXPECT_DOUBLE_EQ(pos[0], 1.0);
    EXPECT_DOUBLE_EQ(pos[1], 2.0);
    EXPECT_DOUBLE_EQ(pos[2], 3.0);

    // Test velocity method
    auto vel = p1.velocity();
    EXPECT_DOUBLE_EQ(vel[0], 0.5);
    EXPECT_DOUBLE_EQ(vel[1], -0.5);
    EXPECT_DOUBLE_EQ(vel[2], 1.0);

    // Test position setting
    p1.set_position({4.0, 5.0, 6.0});
    pos = p1.position();
    EXPECT_DOUBLE_EQ(pos[0], 4.0);
    EXPECT_DOUBLE_EQ(pos[1], 5.0);
    EXPECT_DOUBLE_EQ(pos[2], 6.0);
}

TEST_F(SystemTest, DynamicsUpdate) {
    // Create a valid particle
    Particle p1 = create_test_particle(1, "H1", 0.0, 0.0, 0.0, 0.5);
    p1.set_velocity({1.0, 0.0, 0.0});
    
    // Create a valid residue
    std::vector<Particle> particles = {p1};
    Residue res = create_test_residue(particles, "HOH", 1, 'A');
    
    // Add to system
    system->add_residue(res);
    
    // Update position
    double dt = 0.1;
    system->update_positions(dt);
    
    // Verify position update
    const auto& updated_p1 = system->get_residue(0).atoms[0];
    auto pos = updated_p1.position();
    EXPECT_DOUBLE_EQ(pos[0], dt);
    EXPECT_DOUBLE_EQ(pos[1], 0.0);
    EXPECT_DOUBLE_EQ(pos[2], 0.0);
}

TEST_F(SystemTest, ResidueCenterOfMass) {
    // Create a residue with two particles
    Particle p1 = create_test_particle(1, "H1", 0.0, 0.0, 0.0);
    Particle p2 = create_test_particle(2, "O1", 2.0, 0.0, 0.0);
    Residue res("HOH", 1, 'A');
    res.atoms.push_back(p1);
    res.atoms.push_back(p2);

    // Test center of mass calculation
    auto com = res.center_of_mass();
    EXPECT_DOUBLE_EQ(com[0], 1.0);  // Average x = (0 + 2)/2
    EXPECT_DOUBLE_EQ(com[1], 0.0);
    EXPECT_DOUBLE_EQ(com[2], 0.0);

    // Test atom count
    EXPECT_EQ(res.atom_count(), 2);
}

TEST_F(SystemTest, ParticleVelocityUpdate) {
    // Test velocity updates with different values
    Particle p1 = create_test_particle(1, "H1", 0.0, 0.0, 0.0);
    
    // Test zero velocity
    verify_position_velocity(p1, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0});
    
    // Test setting velocity
    p1.set_velocity({1.0, -1.0, 0.5});
    verify_position_velocity(p1, {0.0, 0.0, 0.0}, {1.0, -1.0, 0.5});
    
    // Test updating position with velocity
    p1.set_position({p1.x + p1.vx, p1.y + p1.vy, p1.z + p1.vz});
    verify_position_velocity(p1, {1.0, -1.0, 0.5}, {1.0, -1.0, 0.5});
}

TEST_F(SystemTest, ResidueCenterOfMassComplex) {
    // Create particles with valid data
    std::vector<Particle> particles = {
        create_test_particle(1, "H1", -2.0, -2.0, -2.0),  // 第一个粒子
        create_test_particle(2, "O1", 2.0, 2.0, 2.0),     // 第二个粒子
        create_test_particle(3, "H2", 4.0, 4.0, 4.0)      // 第三个粒子
    };
    
    Residue res = create_test_residue(particles, "HOH", 1, 'A');
    
    // 计算期望的质心：(-2+2+4)/3 = 4/3 ≈ 1.33
    verify_center_of_mass(res, std::array<double, 3>{4.0/3.0, 4.0/3.0, 4.0/3.0});
    
    // 移动第一个粒子
    res.atoms[0].set_position({-4.0, -4.0, -4.0});
    // 新的期望质心：(-4+2+4)/3 = 2/3 ≈ 0.67
    verify_center_of_mass(res, std::array<double, 3>{2.0/3.0, 2.0/3.0, 2.0/3.0});
}

TEST_F(SystemTest, VelocityIntegration) {
    // Test velocity integration over multiple timesteps
    Particle p1 = create_test_particle(1, "H1", 0.0, 0.0, 0.0);
    p1.set_velocity({1.0, 0.5, 0.0});
    
    Residue res = create_test_residue({p1}, "TEST", 1, 'A');
    system->add_residue(res);
    
    double dt = 0.1;
    int steps = 5;
    
    for (int i = 0; i < steps; ++i) {
        system->update_positions(dt);
        auto& particle = system->get_residue(0).atoms[0];
        auto pos = particle.position();
        
        // Position should increase linearly with time
        EXPECT_DOUBLE_EQ(pos[0], (i + 1) * dt * 1.0);
        EXPECT_DOUBLE_EQ(pos[1], (i + 1) * dt * 0.5);
        EXPECT_DOUBLE_EQ(pos[2], 0.0);
    }
}

TEST_F(SystemTest, ResidueGeometry) {
    // Create a water molecule with realistic geometry
    std::vector<Particle> water = {
        create_test_particle(1, "O", 0.0, 0.0, 0.0),
        create_test_particle(2, "H1", 0.957, 0.0, 0.0),
        create_test_particle(3, "H2", -0.24, 0.927, 0.0)
    };
    
    Residue res = create_test_residue(water, "HOH", 1, 'A');
    system->add_residue(res);
    
    // Verify center of mass
    auto expected_com = std::array<double, 3>{
        (0.0 + 0.957 - 0.24) / 3.0,
        (0.0 + 0.0 + 0.927) / 3.0,
        0.0
    };
    verify_center_of_mass(res, expected_com);
    
    // Test rotation by updating positions
    for (auto& atom : res.atoms) {
        auto pos = atom.position();
        atom.set_position({-pos[1], pos[0], pos[2]});  // 90-degree rotation around z
    }
    
    // Verify rotated center of mass (should be same magnitude, different direction)
    auto com = res.center_of_mass();
    EXPECT_NEAR(std::hypot(com[0], com[1]), 
                std::hypot(expected_com[0], expected_com[1]), 
                1e-10);
}

} // namespace
