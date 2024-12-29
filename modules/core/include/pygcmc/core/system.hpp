// modules/core/include/pygcmc/core/system.hpp

#ifndef PYGCMC_CORE_SYSTEM_HPP
#define PYGCMC_CORE_SYSTEM_HPP

#include "pygcmc/core/io/parser_common.hpp"
#include "pygcmc/core/io/pdb_parser.hpp"
#include "pygcmc/core/io/psf_parser.hpp"
#include "pygcmc/core/io/top_parser.hpp"
#include "pygcmc/core/io/itp_parser.hpp"
#include "pygcmc/core/io/ff_parser.hpp"
#include <vector>
#include <string>
#include <array>
#include <cmath>
#include <stdexcept>

namespace pygcmc {
namespace core {

/**
 * @brief System-related exceptions
 */
class SystemError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/**
 * @brief Represents a particle in the system
 */
struct Particle {
    int serial;              ///< Atom serial number
    std::string name;        ///< Atom name
    std::string residue;     ///< Residue name
    int sequence;            ///< Residue sequence number
    double x, y, z;          ///< Atomic coordinates
    double charge;           ///< Atomic charge
    std::string type;        ///< Atom type
    std::string nameTop;     ///< Topology name
    int typeNum;             ///< Atom type number
    double vx, vy, vz;       ///< Velocities

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

    void set_velocity(const std::array<double, 3>& v) { 
        vx = v[0]; 
        vy = v[1]; 
        vz = v[2]; 
    }

    void set_position(const std::array<double, 3>& pos) { 
        x = pos[0]; 
        y = pos[1]; 
        z = pos[2]; 
    }

    // Method to apply periodic boundary conditions
    void apply_periodic_boundary(double box_size) {
        x = x - box_size * std::floor(x / box_size);
        y = y - box_size * std::floor(y / box_size);
        z = z - box_size * std::floor(z / box_size);
    }
};

/**
 * @brief Represents a residue in the system
 */
struct Residue {
    std::string name;                ///< Residue name
    int sequence_number;             ///< Residue sequence number
    char chain_id;                   ///< Chain identifier
    std::vector<Particle> atoms;      ///< Particles within the residue

    Residue(const std::string& name_ = "", int seq_num_ = 0, char chain_ = ' ')
        : name(name_), sequence_number(seq_num_), chain_id(chain_) {}

    bool is_valid() const {
        if (name.empty()) return false;
        if (sequence_number <= 0) return false;
        if (atoms.empty()) return false;
        for (const auto& atom : atoms) {
            if (!atom.is_valid()) return false;
        }
        return true;
    }

    size_t get_atom_count() const {
        return atoms.size();
    }

    // Compute center of mass (assuming equal mass)
    std::array<double, 3> center_of_mass() const {
        if (atoms.empty()) {
            return {0.0, 0.0, 0.0};
        }
        double sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
        for (const auto& atom : atoms) {
            sum_x += atom.x;
            sum_y += atom.y;
            sum_z += atom.z;
        }
        double n = static_cast<double>(atoms.size());
        return {sum_x / n, sum_y / n, sum_z / n};
    }

    // Get number of atoms
    size_t atom_count() const { 
        return atoms.size(); 
    }
};

/**
 * @brief Manages the molecular system
 */
class System {
public:
    /**
     * @brief Constructs the system with given force field parameters
     * @param epsilon Lennard-Jones epsilon parameter
     * @param sigma Lennard-Jones sigma parameter
     */
    System(double epsilon = 1.0, double sigma = 1.0);
    
    ~System();

    // File loading methods
    void load_pdb(const std::string& filename);
    void load_psf(const std::string& filename);
    void load_top(const std::string& filename);
    void load_itp(const std::string& filename);
    void load_forcefield(const std::string& filename);

    // Residue management
    void add_residue(const Residue& residue);
    void remove_residue(int index);
    size_t get_residue_count() const;
    const Residue& get_residue(size_t index) const;
    Residue& get_residue(size_t index);

    // Energy computation
    double compute_total_energy() const;
    std::pair<double, double> get_system_state() const;

    // Dynamics methods
    void update_positions(double dt);
    void update_velocities(double dt);

    // Boundary conditions
    void set_periodic_boundary(double box_size);
    double apply_pbc(double x) const;

    void add_particle(const Particle& particle) {
        // Create a new residue for the particle
        Residue res(particle.residue, particle.sequence);
        res.atoms.push_back(particle);
        residues_.push_back(res);
    }

    void remove_particle(size_t index) {
        size_t current = 0;
        for (auto it = residues_.begin(); it != residues_.end(); ++it) {
            if (current + it->atoms.size() > index) {
                size_t local_index = index - current;
                it->atoms.erase(it->atoms.begin() + local_index);
                if (it->atoms.empty()) {
                    residues_.erase(it);
                }
                return;
            }
            current += it->atoms.size();
        }
    }

    size_t get_particle_count() const {
        size_t count = 0;
        for (const auto& res : residues_) {
            count += res.atoms.size();
        }
        return count;
    }

    // Particle access methods
    const Particle& get_particle(size_t index) const {
        size_t current = 0;
        for (const auto& res : residues_) {
            if (current + res.atoms.size() > index) {
                return res.atoms[index - current];
            }
            current += res.atoms.size();
        }
        throw std::out_of_range("Particle index out of range");
    }

    Particle& get_particle(size_t index) {
        size_t current = 0;
        for (auto& res : residues_) {
            if (current + res.atoms.size() > index) {
                return res.atoms[index - current];
            }
            current += res.atoms.size();
        }
        throw std::out_of_range("Particle index out of range");
    }

private:
    std::vector<Residue> residues_; ///< List of residues in the system
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
    double compute_residue_energy(const Residue& res1, const Residue& res2) const;
};

} // namespace core
} // namespace pygcmc

#endif // PYGCMC_CORE_SYSTEM_HPP
