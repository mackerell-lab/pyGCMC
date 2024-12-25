// tests/core/test_system.cpp
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include "pygcmc/core/system.hpp"

TEST_CASE("System class basic functionality", "[System]") {
    System system;

    REQUIRE(system.get_particle_count() == 0);

    Particle p1(1.0, 2.0, 3.0, 1, 0.5);
    system.add_particle(p1);
    REQUIRE(system.get_particle_count() == 1);

    Particle p2(4.0, 5.0, 6.0, 2, -0.5);
    system.add_particle(p2);
    REQUIRE(system.get_particle_count() == 2);

    double energy = system.compute_total_energy();
    REQUIRE(energy != Approx(0.0)); // 根据粒子位置和LJ势计算具体值

    system.remove_particle(0);
    REQUIRE(system.get_particle_count() == 1);

    system.remove_particle(10); // 超出范围，应该不改变
    REQUIRE(system.get_particle_count() == 1);
}
