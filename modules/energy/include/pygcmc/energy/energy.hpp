#pragma once

#include <pygcmc/core/project.hpp>
#include <array>

namespace pygcmc {

class Energy {
public:
    virtual ~Energy() = default;
    
    // Set periodic box dimensions
    void set_box(const std::array<double, 3>& box_dimensions) {
        box_ = box_dimensions;
    }
    
    // Get periodic box dimensions
    const std::array<double, 3>& get_box() const {
        return box_;
    }
    
    // Calculate energy for the entire system
    virtual double calculate(const Project& project) = 0;
    
    // Calculate energy for a specific atom with the rest of the system
    virtual double calculate_atom_energy(const Project& project, const ProjectAtom& atom) = 0;
    
protected:
    Energy() = default;
    
    // Apply minimum image convention
    void apply_minimum_image(double& dx, double& dy, double& dz) const {
        dx -= box_[0] * std::round(dx / box_[0]);
        dy -= box_[1] * std::round(dy / box_[1]);
        dz -= box_[2] * std::round(dz / box_[2]);
    }
    
private:
    std::array<double, 3> box_{1.0e10, 1.0e10, 1.0e10}; // Default to very large box (effectively no PBC)
};

} // namespace pygcmc 