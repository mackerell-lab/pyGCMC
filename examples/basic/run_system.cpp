// examples/basic/run_system.cpp
#include <iostream>
#include <iomanip>
#include "pygcmc/core/system.hpp"

using pygcmc::core::System;
using pygcmc::core::Particle;
using pygcmc::core::Residue;

void print_particle_info(const Particle& p) {
    auto pos = p.position();
    std::cout << "Particle " << p.serial << " (" << p.name << "):\n"
              << "  Position: (" << std::fixed << std::setprecision(3)
              << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n"
              << "  Charge: " << p.charge << "\n"
              << "  Type: " << p.type << "\n";
}

void print_system_state(const System& sys) {
    auto [kinetic, potential] = sys.get_system_state();
    std::cout << "\nSystem State:\n"
              << "  Residue count: " << sys.get_residue_count() << "\n"
              << "  Kinetic energy: " << std::scientific << kinetic << "\n"
              << "  Potential energy: " << potential << "\n"
              << "  Total energy: " << kinetic + potential << "\n"
              << std::defaultfloat;

    // Print center of mass for each residue
    for (size_t i = 0; i < sys.get_residue_count(); ++i) {
        const auto& res = sys.get_residue(i);
        auto com = res.center_of_mass();
        std::cout << "  Residue " << i << " center: ("
                  << std::fixed << std::setprecision(3)
                  << com[0] << ", " << com[1] << ", " << com[2] << ")\n";
    }
}

int main() {
    try {
        // Initialize system with custom LJ parameters
        System sys(1.0, 3.355);  // epsilon = 1.0, sigma = 3.355

        // Create water molecule particles
        Particle h1(1, "H1", "HOH", 1, 0.0, 0.0, 0.0, 0.5, "H", "H");
        Particle o1(2, "O1", "HOH", 1, 1.0, 0.0, 0.0, -1.0, "O", "O");
        Particle h2(3, "H2", "HOH", 1, 1.0, 1.0, 0.0, 0.5, "H", "H");

        // Create a residue and add particles to it
        Residue water("HOH", 1, 'A');
        water.atoms.push_back(h1);
        water.atoms.push_back(o1);
        water.atoms.push_back(h2);

        // Add residue to system
        sys.add_residue(water);

        // Print initial state
        std::cout << "Initial configuration:\n";
        for (size_t i = 0; i < sys.get_residue_count(); ++i) {
            const Residue& res = sys.get_residue(i);
            std::cout << "Residue " << res.sequence_number << " (" << res.name 
                      << ") with " << res.atom_count() << " atoms:\n";
            for (const auto& atom : res.atoms) {
                print_particle_info(atom);
            }
        }
        print_system_state(sys);

        // Set up periodic boundary conditions
        sys.set_periodic_boundary(10.0);

        // Add some initial velocities
        auto& res0 = sys.get_residue(0);
        if (!res0.atoms.empty()) {
            auto& p1 = res0.atoms[0];
            p1.set_velocity({1.0, 0.0, 0.0});  // Using new set_velocity method
        }

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
