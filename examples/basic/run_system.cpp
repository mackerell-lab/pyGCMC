// examples/basic/run_system.cpp
#include <iostream>
#include <iomanip>
#include "pygcmc/core/system.hpp"

using pygcmc::core::System;
using pygcmc::core::Particle;

void print_particle_info(const Particle& p) {
    std::cout << "Particle " << p.serial << " (" << p.name << "):\n"
              << "  Position: (" << std::fixed << std::setprecision(3)
              << p.x << ", " << p.y << ", " << p.z << ")\n"
              << "  Charge: " << p.charge << "\n"
              << "  Type: " << p.type << "\n";
}

void print_system_state(const System& sys) {
    auto [kinetic, potential] = sys.get_system_state();
    std::cout << "\nSystem State:\n"
              << "  Particle count: " << sys.get_particle_count() << "\n"
              << "  Kinetic energy: " << std::scientific << kinetic << "\n"
              << "  Potential energy: " << potential << "\n"
              << "  Total energy: " << kinetic + potential << "\n"
              << std::defaultfloat;
}

int main() {
    try {
        // Initialize system with custom LJ parameters
        System sys(1.0, 3.355);  // epsilon = 1.0, sigma = 3.355

        // Create water molecule particles
        Particle h1(1, "H1", "HOH", 1, 0.0, 0.0, 0.0, 0.5, 1, "H");
        Particle o1(2, "O1", "HOH", 1, 1.0, 0.0, 0.0, -1.0, 2, "O");
        Particle h2(3, "H2", "HOH", 1, 1.0, 1.0, 0.0, 0.5, 1, "H");

        // Add particles to system
        sys.add_particle(h1);
        sys.add_particle(o1);
        sys.add_particle(h2);

        // Print initial state
        std::cout << "Initial configuration:\n";
        for (size_t i = 0; i < sys.get_particle_count(); ++i) {
            print_particle_info(sys.get_particle(i));
        }
        print_system_state(sys);

        // Set up periodic boundary conditions
        sys.set_periodic_boundary(10.0);

        // Add some initial velocities
        auto& p1 = sys.get_particle(0);
        p1.set_velocity(1.0, 0.0, 0.0);

        // Run a few dynamics steps
        double dt = 0.001;
        for (int step = 0; step < 10; ++step) {
            sys.update_positions(dt);
            sys.update_velocities(dt);
            
            if (step % 5 == 0) {
                std::cout << "\nStep " << step << ":\n";
                print_system_state(sys);
            }
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
