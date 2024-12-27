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
        return Particle(serial, name, "TEST", 1, x, y, z, charge);
    }

    std::unique_ptr<System> system;
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

TEST_F(SystemTest, DynamicsUpdate) {
    auto p1 = create_test_particle(1, "H1", 0.0, 0.0, 0.0);
    p1.vx = 1.0;
    system->add_particle(p1);
    
    // Update position
    double dt = 0.1;
    system->update_positions(dt);
    
    const auto& updated_particle = system->get_particle(0);
    EXPECT_DOUBLE_EQ(updated_particle.x, dt);
    EXPECT_DOUBLE_EQ(updated_particle.y, 0.0);
    EXPECT_DOUBLE_EQ(updated_particle.z, 0.0);
}

} // namespace
