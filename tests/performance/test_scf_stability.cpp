// Test SCF stability with detailed output

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

const double ONE_4PI_EPS0 = 138.935456;
const double ANGSTROM_TO_NM = 0.1;
const double KCAL_TO_KJ = 4.184;

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(double s) const { return Vec3(x * s, y * s, z * s); }
    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    double norm2() const { return x*x + y*y + z*z; }
};

void testSCFWithPreconditioning() {
    std::cout << "=== Testing SCF with Preconditioning ===\n\n";
    
    // Parameters
    const double qO = 1.71636;
    const double qD = -1.71636;
    const double qH = 0.55733;
    const double qM = -1.11466;
    const double k_drude = 1000.0 * KCAL_TO_KJ / (ANGSTROM_TO_NM * ANGSTROM_TO_NM);
    
    // Create three water molecules in a line
    std::vector<Vec3> pos;
    std::vector<double> charges;
    std::vector<std::string> names;
    std::vector<int> molecule;
    std::vector<int> drude_indices;
    std::vector<int> parent_indices;
    
    double spacing = 0.3;  // nm between molecules
    
    for (int mol = 0; mol < 3; mol++) {
        double x_offset = mol * spacing;
        int idx_start = pos.size();
        
        // Oxygen
        pos.push_back(Vec3(x_offset, 0, 0));
        charges.push_back(qO);
        names.push_back("O" + std::to_string(mol+1));
        molecule.push_back(mol);
        
        // Drude - start VERY close to parent
        pos.push_back(Vec3(x_offset, 0, 1e-8));  // 1e-8 nm = 0.0001 Å
        charges.push_back(qD);
        names.push_back("D" + std::to_string(mol+1));
        molecule.push_back(mol);
        drude_indices.push_back(idx_start + 1);
        parent_indices.push_back(idx_start);
        
        // Simplified - just add one H for clarity
        pos.push_back(Vec3(x_offset + 0.09572, 0, 0));
        charges.push_back(qH);
        names.push_back("H" + std::to_string(mol+1));
        molecule.push_back(mol);
    }
    
    std::cout << "Initial configuration:\n";
    for (size_t i = 0; i < drude_indices.size(); i++) {
        Vec3 dr = pos[drude_indices[i]] - pos[parent_indices[i]];
        std::cout << "  " << names[drude_indices[i]] << " displacement: " 
                  << dr.norm() << " nm\n";
    }
    
    // Step 1: Precondition - optimize each Drude in isolation
    std::cout << "\nStep 1: Preconditioning (isolated optimization)...\n";
    
    for (size_t d = 0; d < drude_indices.size(); d++) {
        int di = drude_indices[d];
        int pi = parent_indices[d];
        
        // For isolated Drude, equilibrium is at parent position
        pos[di] = pos[pi];
        std::cout << "  " << names[di] << " moved to parent position\n";
    }
    
    // Step 2: Full SCF with all interactions
    std::cout << "\nStep 2: Full SCF optimization...\n";
    
    double tolerance = 0.1;  // kJ/mol/nm
    int max_iter = 100;
    
    for (int iter = 0; iter < max_iter; iter++) {
        double max_force = 0.0;
        std::vector<Vec3> forces(drude_indices.size());
        
        // Calculate all forces first
        for (size_t d = 0; d < drude_indices.size(); d++) {
            int di = drude_indices[d];
            int pi = parent_indices[d];
            
            Vec3 force(0, 0, 0);
            
            // Harmonic force
            Vec3 dr = pos[di] - pos[pi];
            force = force - (dr * k_drude);
            
            // Coulomb forces from other atoms (not in same molecule)
            for (size_t j = 0; j < pos.size(); j++) {
                if (molecule[j] != molecule[di]) {
                    Vec3 rij = pos[j] - pos[di];
                    double r = rij.norm();
                    
                    if (r > 1e-10) {
                        double f_mag = ONE_4PI_EPS0 * charges[di] * charges[j] / (r * r);
                        force = force - (rij * (f_mag / r));
                    }
                }
            }
            
            forces[d] = force;
            max_force = std::max(max_force, force.norm());
        }
        
        if (iter % 10 == 0 || max_force < tolerance) {
            std::cout << "  Iter " << iter << ": max force = " << max_force << " kJ/mol/nm\n";
            if (iter == 0) {
                for (size_t d = 0; d < drude_indices.size(); d++) {
                    std::cout << "    " << names[drude_indices[d]] << " force: " 
                              << forces[d].norm() << " kJ/mol/nm\n";
                }
            }
        }
        
        if (max_force < tolerance) {
            std::cout << "  Converged after " << iter << " iterations\n";
            break;
        }
        
        // Update all Drude positions
        for (size_t d = 0; d < drude_indices.size(); d++) {
            if (forces[d].norm() > tolerance * 0.01) {
                // Very conservative update
                double damping = 0.2;
                double force_mag = forces[d].norm();
                
                // Reduce damping for large forces
                if (force_mag > 1000.0) {
                    damping = 0.001;
                } else if (force_mag > 100.0) {
                    damping = 0.01;
                } else if (force_mag > 10.0) {
                    damping = 0.1;
                }
                
                Vec3 delta = forces[d] * (damping / k_drude);
                
                // Hard limit on displacement
                double max_disp = 0.0001;  // 0.0001 nm = 0.001 Å
                double delta_mag = delta.norm();
                if (delta_mag > max_disp) {
                    delta = delta * (max_disp / delta_mag);
                }
                
                pos[drude_indices[d]] = pos[drude_indices[d]] + delta;
            }
        }
    }
    
    // Final results
    std::cout << "\nFinal Drude displacements:\n";
    double total_dipole = 0.0;
    for (size_t d = 0; d < drude_indices.size(); d++) {
        Vec3 dr = pos[drude_indices[d]] - pos[parent_indices[d]];
        double disp = dr.norm();
        std::cout << "  " << names[drude_indices[d]] << ": " << disp << " nm";
        if (disp > 0.01) {
            std::cout << " (WARNING: Large displacement!)";
        }
        std::cout << "\n";
        total_dipole += disp * std::abs(qD);
    }
    
    std::cout << "\nAverage induced dipole: " << total_dipole/drude_indices.size() 
              << " e·nm\n";
    
    // Calculate final energy
    double e_harmonic = 0.0;
    for (size_t d = 0; d < drude_indices.size(); d++) {
        Vec3 dr = pos[drude_indices[d]] - pos[parent_indices[d]];
        e_harmonic += 0.5 * k_drude * dr.norm2();
    }
    std::cout << "\nTotal harmonic energy: " << e_harmonic << " kJ/mol\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(8);
    testSCFWithPreconditioning();
    return 0;
}