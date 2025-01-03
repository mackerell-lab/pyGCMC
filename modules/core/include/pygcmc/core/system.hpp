// modules/core/include/pygcmc/core/system.hpp

#ifndef PYGCMC_CORE_SYSTEM_HPP
#define PYGCMC_CORE_SYSTEM_HPP

#include <vector>
#include <array>
#include <memory>
#include <string>
#include <cmath>
#include <algorithm>
#include "pygcmc/core/force.hpp"
#include "pygcmc/core/structure.hpp"

namespace pygcmc {
namespace core {

/**
 * @brief Represents a particle in the system
 */
struct Particle {
    std::array<double, 3> position;
    std::array<double, 3> velocity;
    double charge;
    double mass;
    bool is_virtual;

    Particle(const std::array<double, 3>& pos = {0, 0, 0},
            const std::array<double, 3>& vel = {0, 0, 0},
            double q = 0.0,
            double m = 1.0)
        : position(pos), velocity(vel), charge(q), mass(m), is_virtual(false) {}

    bool is_valid() const {
        return std::isfinite(charge) && std::isfinite(mass) && mass >= 0.0 &&
               std::all_of(position.begin(), position.end(), [](double x) { return std::isfinite(x); }) &&
               std::all_of(velocity.begin(), velocity.end(), [](double x) { return std::isfinite(x); });
    }
};

/**
 * @brief Represents a constraint between two particles
 */
struct Constraint {
    size_t particle1;
    size_t particle2;
    double distance;

    Constraint(size_t p1, size_t p2, double d)
        : particle1(p1), particle2(p2), distance(d) {}
};

/**
 * @brief Represents a residue in the system
 */
struct Residue {
    std::string name;
    std::vector<Particle> particles;

    Residue(const std::string& n = "") : name(n) {}

    size_t atom_count() const { return particles.size(); }
    
    std::array<double, 3> center_of_mass() const {
        if (particles.empty()) {
            return {0.0, 0.0, 0.0};
        }
        
        double total_mass = 0.0;
        std::array<double, 3> com = {0.0, 0.0, 0.0};
        
        for (const auto& particle : particles) {
            if (!particle.is_virtual) {
                total_mass += particle.mass;
                for (size_t i = 0; i < 3; ++i) {
                    com[i] += particle.position[i] * particle.mass;
                }
            }
        }
        
        if (total_mass > 0.0) {
            for (size_t i = 0; i < 3; ++i) {
                com[i] /= total_mass;
            }
        }
        
        return com;
    }
};

/**
 * @brief Main class representing a molecular system
 */
class System {
public:
    // Default constructor
    System() = default;
    
    // Constructor with PDB and PSF files
    System(const std::string& pdb, const std::string& psf) {
        if (!pdb.empty() && !psf.empty()) {
            load_structure_psf(pdb, psf);
        }
    }
    
    // Static factory method for TOP files
    static System from_top(const std::string& pdb, const std::string& top) {
        System system;
        system.load_structure_top(pdb, top);
        return system;
    }
    
    ~System() = default;

    // Structure loading methods
    void load_structure_psf(const std::string& pdb, const std::string& psf);
    void load_structure_top(const std::string& pdb, const std::string& top);
    void load_structure(const Structure& structure);

    // Residue management
    size_t add_residue(const std::string& name);
    void remove_residue(size_t index);
    size_t get_residue_count() const;
    const Residue& get_residue(size_t index) const;
    Residue& get_residue(size_t index);

    // Particle management
    size_t add_particle(size_t residue_index, const Particle& particle);
    void remove_particle(size_t residue_index, size_t particle_index);
    size_t get_particle_count(size_t residue_index) const;
    const Particle& get_particle(size_t residue_index, size_t particle_index) const;
    Particle& get_particle(size_t residue_index, size_t particle_index);

    // Mass management
    double get_particle_mass(size_t residue_index, size_t particle_index) const;
    void set_particle_mass(size_t residue_index, size_t particle_index, double mass);

    // Virtual site management
    void set_virtual_site(size_t residue_index, size_t particle_index, bool is_virtual);
    bool is_virtual_site(size_t residue_index, size_t particle_index) const;

    // Constraint management
    size_t add_constraint(size_t residue1, size_t particle1, size_t residue2, size_t particle2, double distance);
    void remove_constraint(size_t index);
    size_t get_constraint_count() const;
    const Constraint& get_constraint(size_t index) const;
    void get_constraint_parameters(size_t index, size_t& residue1, size_t& particle1, 
                                 size_t& residue2, size_t& particle2, double& distance) const;
    void set_constraint_parameters(size_t index, size_t residue1, size_t particle1,
                                 size_t residue2, size_t particle2, double distance);

    // Force management
    void add_force(std::shared_ptr<Force> force);
    void remove_force(size_t index);
    size_t get_force_count() const;
    std::shared_ptr<Force> get_force(size_t index) const;

    // Periodic boundary conditions
    void set_periodic_box_vectors(const std::array<double, 3>& a,
                                const std::array<double, 3>& b,
                                const std::array<double, 3>& c);
    void get_periodic_box_vectors(std::array<double, 3>& a,
                                std::array<double, 3>& b,
                                std::array<double, 3>& c) const;
    bool uses_periodic_boundary_conditions() const;

    // Distance computation
    std::array<double, 3> compute_distance(const Particle& p1, const Particle& p2) const;

    // System state
    double compute_energy() const;
    void update_positions(double dt);
    void update_velocities(double dt);

private:
    std::vector<Residue> residues_;
    std::vector<Constraint> constraints_;
    std::vector<std::shared_ptr<Force>> forces_;
    std::array<std::array<double, 3>, 3> box_vectors_;
    bool has_periodic_boundary_ = false;

    void validate_residue_index(size_t index) const;
    void validate_particle_index(size_t residue_index, size_t particle_index) const;
    void validate_constraint_index(size_t index) const;
    void validate_force_index(size_t index) const;
    void validate_box_vectors(const std::array<double, 3>& a,
                            const std::array<double, 3>& b,
                            const std::array<double, 3>& c) const;
};

} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_SYSTEM_HPP
