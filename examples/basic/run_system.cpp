// examples/basic/run_system.cpp
#include <iostream>
#include <array>
#include <string>
#include "pygcmc/core/system.hpp"

using namespace pygcmc::core;

void print_particle_info(const Particle& p) {
    std::cout << "Particle:\n"
              << "  Position: [" << p.position[0] << ", " << p.position[1] << ", " << p.position[2] << "]\n"
              << "  Velocity: [" << p.velocity[0] << ", " << p.velocity[1] << ", " << p.velocity[2] << "]\n"
              << "  Charge: " << p.charge << "\n"
              << "  Mass: " << p.mass << "\n"
              << "  Virtual: " << (p.is_virtual ? "yes" : "no") << "\n";
}

void print_system_state(const System& sys) {
    double energy = sys.compute_energy();
    std::cout << "System energy: " << energy << "\n";
}

int main() {
    try {
        // Create a new system
        System sys;

        // Create particles for a water molecule
        Particle h1({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.5, 1.008);  // H1
        Particle o1({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, -1.0, 15.999); // O
        Particle h2({1.0, 1.0, 0.0}, {0.0, 0.0, 0.0}, 0.5, 1.008);  // H2

        // Create a residue for water
        size_t water_idx = sys.add_residue("HOH");
        Residue& water = sys.get_residue(water_idx);

        // Add particles to the water residue
        sys.add_particle(water_idx, h1);
        sys.add_particle(water_idx, o1);
        sys.add_particle(water_idx, h2);

        // Print initial system state
        std::cout << "Initial system state:\n";
        for (size_t i = 0; i < sys.get_residue_count(); ++i) {
            const Residue& res = sys.get_residue(i);
            std::cout << "Residue " << i << " (" << res.name << "):\n";
            
            for (size_t j = 0; j < sys.get_particle_count(i); ++j) {
                const Particle& particle = sys.get_particle(i, j);
                print_particle_info(particle);
            }
        }

        // Set up periodic boundary conditions
        std::array<double, 3> a = {10.0, 0.0, 0.0};
        std::array<double, 3> b = {0.0, 10.0, 0.0};
        std::array<double, 3> c = {0.0, 0.0, 10.0};
        sys.set_periodic_box_vectors(a, b, c);

        // Add constraints between O-H bonds
        sys.add_constraint(water_idx, 1, water_idx, 0, 0.9572); // O-H1 bond
        sys.add_constraint(water_idx, 1, water_idx, 2, 0.9572); // O-H2 bond

        // Print final system state
        std::cout << "\nFinal system state:\n";
        print_system_state(sys);

        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
