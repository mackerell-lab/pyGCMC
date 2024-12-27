// modules/core/include/pygcmc/core/system.hpp

#ifndef PYGCMC_CORE_SYSTEM_HPP
#define PYGCMC_CORE_SYSTEM_HPP

#include "pygcmc/core/io/parser_common.hpp"
#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <stdexcept>
#include <memory>

namespace pygcmc {
namespace core {

/**
 * @brief System-related exceptions
 */
class SystemError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

/**
 * @brief Represents a particle in the system
 */
struct Particle {
    int serial;
    std::string name;
    std::string residue;
    int sequence;
    double x, y, z;
    double charge;
    std::string type;
    std::string nameTop;
    int typeNum;
    double vx, vy, vz;

    Particle(int serial_ = 0, const std::string& name_ = "", 
             const std::string& residue_ = "", int sequence_ = 0,
             double x_ = 0.0, double y_ = 0.0, double z_ = 0.0,
             double charge_ = 0.0, const std::string& type_ = "", 
             const std::string& nameTop_ = "")
        : serial(serial_), name(name_), residue(residue_), 
          sequence(sequence_), x(x_), y(y_), z(z_), 
          charge(charge_), type(type_), nameTop(nameTop_),
          typeNum(0), vx(0.0), vy(0.0), vz(0.0) {}

    bool is_valid() const {
        return serial > 0 && 
               !name.empty() && 
               !residue.empty() && 
               sequence > 0 &&
               std::isfinite(x) && std::isfinite(y) && std::isfinite(z) &&
               std::isfinite(charge) &&
               std::isfinite(vx) && std::isfinite(vy) && std::isfinite(vz) &&
               !type.empty();
    }

    std::array<double, 3> position() const { 
        return {x, y, z}; 
    }

    std::array<double, 3> velocity() const { 
        return {vx, vy, vz}; 
    }
    
    void set_position(const std::array<double, 3>& pos) {
        x = pos[0]; y = pos[1]; z = pos[2];
    }
    
    void set_velocity(const std::array<double, 3>& vel) {
        vx = vel[0]; vy = vel[1]; vz = vel[2];
    }

    double kinetic_energy() const {
        return 0.5 * (vx * vx + vy * vy + vz * vz);
    }

    void apply_periodic_boundary(double box_size) {
        x -= box_size * std::floor(x / box_size);
        y -= box_size * std::floor(y / box_size);
        z -= box_size * std::floor(z / box_size);
    }
};

class System {
public:
    System(double epsilon = 1.0, double sigma = 1.0);
    ~System();

    // File loading methods
    void load_pdb(const std::string& filename);
    void load_psf(const std::string& filename);
    void load_top(const std::string& filename);
    void load_itp(const std::string& filename);
    void load_forcefield(const std::string& filename);

    // Particle management
    void add_particle(const Particle& particle);
    void remove_particle(int index);
    size_t get_particle_count() const;
    const Particle& get_particle(size_t index) const;
    Particle& get_particle(size_t index);

    // Energy computation
    double compute_total_energy() const;
    std::pair<double, double> get_system_state() const;

    // Dynamics methods
    void update_positions(double dt);
    void update_velocities(double dt);

    // Boundary conditions
    void set_periodic_boundary(double box_size);
    double apply_pbc(double x) const;

private:
    std::vector<Particle> particles_;
    double epsilon_;
    double sigma_;
    io::NBMap nb_dict_;
    io::NBFixMap nbfix_dict_;
    double box_size_ = 0.0;
    bool use_periodic_ = false;

    std::array<double, 3> compute_distance(const Particle& p1, 
                                         const Particle& p2) const;
    double compute_pair_energy(const Particle& p1, 
                             const Particle& p2) const;
};

} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_SYSTEM_HPP
