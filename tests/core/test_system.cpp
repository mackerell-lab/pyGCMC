// tests/core/test_system.cpp
#include <gtest/gtest.h>
#include "pygcmc/core/system.hpp"
#include <memory>
#include <cmath>

using namespace pygcmc::core;

class SystemTest : public ::testing::Test {
protected:
    void SetUp() override {
        system = std::make_unique<System>();
    }

    std::unique_ptr<System> system;
};

// Test residue management
TEST_F(SystemTest, ResidueManagement) {
    size_t res_idx = system->add_residue("ALA");
    EXPECT_EQ(system->get_residue_count(), 1);
    EXPECT_EQ(system->get_residue(res_idx).name, "ALA");
    
    system->remove_residue(res_idx);
    EXPECT_EQ(system->get_residue_count(), 0);
    
    EXPECT_THROW(system->get_residue(0), std::out_of_range);
}

// Test particle management
TEST_F(SystemTest, ParticleManagement) {
    size_t res_idx = system->add_residue("ALA");
    Particle p({1.0, 2.0, 3.0}, {0.1, 0.2, 0.3}, 1.0, 2.0);
    
    size_t p_idx = system->add_particle(res_idx, p);
    EXPECT_EQ(system->get_particle_count(res_idx), 1);
    
    const auto& added_p = system->get_particle(res_idx, p_idx);
    EXPECT_EQ(added_p.position[0], 1.0);
    EXPECT_EQ(added_p.velocity[1], 0.2);
    EXPECT_EQ(added_p.charge, 1.0);
    EXPECT_EQ(added_p.mass, 2.0);
    
    system->remove_particle(res_idx, p_idx);
    EXPECT_EQ(system->get_particle_count(res_idx), 0);
}

// Test mass management
TEST_F(SystemTest, MassManagement) {
    size_t res_idx = system->add_residue("ALA");
    Particle p({0, 0, 0}, {0, 0, 0}, 0.0, 1.0);
    size_t p_idx = system->add_particle(res_idx, p);
    
    EXPECT_EQ(system->get_particle_mass(res_idx, p_idx), 1.0);
    
    system->set_particle_mass(res_idx, p_idx, 2.0);
    EXPECT_EQ(system->get_particle_mass(res_idx, p_idx), 2.0);
    
    // Setting mass to zero should make it a virtual site
    system->set_particle_mass(res_idx, p_idx, 0.0);
    EXPECT_TRUE(system->is_virtual_site(res_idx, p_idx));
    
    EXPECT_THROW(system->set_particle_mass(res_idx, p_idx, -1.0), std::invalid_argument);
}

// Test virtual site management
TEST_F(SystemTest, VirtualSiteManagement) {
    size_t res_idx = system->add_residue("ALA");
    Particle p({0, 0, 0}, {0, 0, 0}, 0.0, 1.0);
    size_t p_idx = system->add_particle(res_idx, p);
    
    EXPECT_FALSE(system->is_virtual_site(res_idx, p_idx));
    
    system->set_virtual_site(res_idx, p_idx, true);
    EXPECT_TRUE(system->is_virtual_site(res_idx, p_idx));
    EXPECT_EQ(system->get_particle_mass(res_idx, p_idx), 0.0);
    
    system->set_virtual_site(res_idx, p_idx, false);
    EXPECT_FALSE(system->is_virtual_site(res_idx, p_idx));
}

// Test constraint management
TEST_F(SystemTest, ConstraintManagement) {
    size_t res1_idx = system->add_residue("ALA");
    size_t res2_idx = system->add_residue("GLY");
    
    Particle p1({0, 0, 0}, {0, 0, 0}, 0.0, 1.0);
    Particle p2({1, 0, 0}, {0, 0, 0}, 0.0, 1.0);
    
    size_t p1_idx = system->add_particle(res1_idx, p1);
    size_t p2_idx = system->add_particle(res2_idx, p2);
    
    size_t c_idx = system->add_constraint(res1_idx, p1_idx, res2_idx, p2_idx, 1.0);
    EXPECT_EQ(system->get_constraint_count(), 1);
    
    size_t r1, r2, p1_out, p2_out;
    double dist;
    system->get_constraint_parameters(c_idx, r1, p1_out, r2, p2_out, dist);
    EXPECT_EQ(r1, res1_idx);
    EXPECT_EQ(r2, res2_idx);
    EXPECT_EQ(p1_out, p1_idx);
    EXPECT_EQ(p2_out, p2_idx);
    EXPECT_EQ(dist, 1.0);
    
    system->set_constraint_parameters(c_idx, res1_idx, p1_idx, res2_idx, p2_idx, 2.0);
    system->get_constraint_parameters(c_idx, r1, p1_out, r2, p2_out, dist);
    EXPECT_EQ(dist, 2.0);
    
    // Test invalid operations
    EXPECT_THROW(system->add_constraint(res1_idx, p1_idx, res2_idx, p2_idx, -1.0),
                 std::invalid_argument);  // negative distance
    
    // Make p1 a virtual site
    system->set_virtual_site(res1_idx, p1_idx, true);
    EXPECT_THROW(system->add_constraint(res1_idx, p1_idx, res2_idx, p2_idx, 1.0),
                 std::invalid_argument);  // virtual site constraint
}

// Test periodic boundary conditions
TEST_F(SystemTest, PeriodicBoundaryConditions) {
    std::array<double, 3> a = {2.0, 0.0, 0.0};
    std::array<double, 3> b = {0.0, 2.0, 0.0};
    std::array<double, 3> c = {0.0, 0.0, 2.0};
    
    system->set_periodic_box_vectors(a, b, c);
    
    std::array<double, 3> a_out, b_out, c_out;
    system->get_periodic_box_vectors(a_out, b_out, c_out);
    
    for (size_t i = 0; i < 3; ++i) {
        EXPECT_EQ(a_out[i], a[i]);
        EXPECT_EQ(b_out[i], b[i]);
        EXPECT_EQ(c_out[i], c[i]);
    }
    
    // Test invalid box vectors
    std::array<double, 3> invalid_a = {0.0, 0.0, 0.0};  // zero vector
    EXPECT_THROW(system->set_periodic_box_vectors(invalid_a, b, c),
                 std::invalid_argument);
}

// Mock Force class for testing
class MockForce : public Force {
public:
    double calculate_forces(const System& system,
                          std::vector<std::vector<std::array<double, 3>>>& forces) const override {
        return 0.0;
    }
    
    bool uses_periodic_boundary_conditions() const override {
        return true;
    }
};

// Test force management
TEST_F(SystemTest, ForceManagement) {
    auto force = std::make_shared<MockForce>();
    system->add_force(force);
    EXPECT_EQ(system->get_force_count(), 1);
    
    auto retrieved_force = system->get_force(0);
    EXPECT_EQ(retrieved_force, force);
    
    system->remove_force(0);
    EXPECT_EQ(system->get_force_count(), 0);
    
    EXPECT_THROW(system->get_force(0), std::out_of_range);
    EXPECT_THROW(system->add_force(nullptr), std::invalid_argument);
}

// Test system dynamics
TEST_F(SystemTest, SystemDynamics) {
    size_t res_idx = system->add_residue("ALA");
    Particle p({0, 0, 0}, {1, 1, 1}, 0.0, 1.0);
    size_t p_idx = system->add_particle(res_idx, p);
    
    system->update_positions(0.1);
    const auto& updated_p = system->get_particle(res_idx, p_idx);
    
    for (size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(updated_p.position[i], 0.1, 1e-10);
    }
    
    // Virtual sites should not move
    system->set_virtual_site(res_idx, p_idx, true);
    system->update_positions(0.1);
    const auto& virtual_p = system->get_particle(res_idx, p_idx);
    
    for (size_t i = 0; i < 3; ++i) {
        EXPECT_NEAR(virtual_p.position[i], 0.1, 1e-10);  // position unchanged
    }
}
