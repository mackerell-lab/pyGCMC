// Debug version to find why SWM4-NDP energy is so high

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

void debugSingleWater() {
    std::cout << "=== Debugging Single SWM4-NDP Water Molecule ===\n\n";
    
    // SWM4-NDP parameters
    const double qO = 1.71636;
    const double qD = -1.71636;
    const double qH = 0.55733;
    const double qM = -1.11466;
    const double k_drude = 1000.0 * KCAL_TO_KJ / (ANGSTROM_TO_NM * ANGSTROM_TO_NM);
    
    std::cout << "Parameters:\n";
    std::cout << "  k_drude = " << k_drude << " kJ/mol/nm²\n";
    std::cout << "  Charges: O=" << qO << ", D=" << qD << ", H=" << qH << ", M=" << qM << "\n";
    std::cout << "  Sum of charges = " << (qO + qD + 2*qH + qM) << " (should be 0)\n\n";
    
    // Create a single water molecule
    std::vector<Vec3> pos;
    std::vector<double> charges;
    std::vector<std::string> names;
    
    // Oxygen at origin
    pos.push_back(Vec3(0, 0, 0));
    charges.push_back(qO);
    names.push_back("O");
    
    // Drude - try different initial displacements
    double drude_disp[] = {0.0001, 0.001, 0.01, 0.1};  // nm
    
    for (int test = 0; test < 4; test++) {
        std::cout << "\nTest " << test+1 << ": Drude displacement = " << drude_disp[test] << " nm\n";
        
        // Reset positions for this test
        pos.clear();
        charges.clear();
        names.clear();
        
        // Oxygen
        pos.push_back(Vec3(0, 0, 0));
        charges.push_back(qO);
        names.push_back("O");
        
        // Drude
        pos.push_back(Vec3(0, 0, drude_disp[test]));
        charges.push_back(qD);
        names.push_back("D");
        
        // Calculate just the harmonic energy
        Vec3 dr = pos[1] - pos[0];
        double r = dr.norm();
        double harmonic_energy = 0.5 * k_drude * r * r;
        
        std::cout << "  Drude-parent distance: " << r << " nm\n";
        std::cout << "  Harmonic energy: " << harmonic_energy << " kJ/mol\n";
        
        // Calculate force on Drude
        Vec3 harmonic_force = dr * (-k_drude);
        std::cout << "  Harmonic force magnitude: " << harmonic_force.norm() << " kJ/mol/nm\n";
        
        // Expected equilibrium: F_harmonic + F_coulomb = 0
        // F_coulomb = k*qO*qD/r² (attractive since opposite charges)
        // At equilibrium: k_drude * r = k_coulomb * qO * |qD| / r²
        // r³ = k_coulomb * qO * |qD| / k_drude
        
        double k_coulomb = ONE_4PI_EPS0;
        double r_eq_cubed = k_coulomb * qO * std::abs(qD) / k_drude;
        double r_eq = std::pow(r_eq_cubed, 1.0/3.0);
        
        std::cout << "  Expected equilibrium distance: " << r_eq << " nm\n";
        std::cout << "  Expected equilibrium distance: " << r_eq/ANGSTROM_TO_NM << " Å\n";
    }
    
    // Now test SCF convergence
    std::cout << "\n=== Testing SCF Convergence ===\n";
    
    // Start with Drude at 0.01 nm displacement
    Vec3 o_pos(0, 0, 0);
    Vec3 d_pos(0, 0, 0.01);
    
    double tolerance = 1.0;  // kJ/mol/nm
    int max_iter = 50;
    
    for (int iter = 0; iter < max_iter; iter++) {
        // Calculate forces
        Vec3 dr = d_pos - o_pos;
        Vec3 f_harmonic = dr * (-k_drude);
        
        // For isolated molecule, only harmonic force
        Vec3 f_total = f_harmonic;
        double f_mag = f_total.norm();
        
        if (iter % 10 == 0) {
            std::cout << "Iter " << iter << ": |F| = " << f_mag << " kJ/mol/nm, r = " << dr.norm() << " nm\n";
        }
        
        if (f_mag < tolerance) {
            std::cout << "Converged after " << iter << " iterations\n";
            break;
        }
        
        // Update position
        double damping = (f_mag > 10.0 * tolerance) ? 0.1 : 0.5;
        Vec3 delta = f_total * (damping / k_drude);
        d_pos = d_pos + delta;
    }
    
    std::cout << "\nFinal Drude position: (" << d_pos.x << ", " << d_pos.y << ", " << d_pos.z << ")\n";
    std::cout << "Final distance: " << (d_pos - o_pos).norm() << " nm\n";
    
    // Test with external field
    std::cout << "\n=== Testing with External Electric Field ===\n";
    
    // Reset Drude position
    d_pos = Vec3(0, 0, 0.001);
    
    // Add a point charge at 1 nm distance
    Vec3 ext_pos(1.0, 0, 0);
    double ext_charge = 1.0;
    
    std::cout << "External charge " << ext_charge << " at (" << ext_pos.x << ", " << ext_pos.y << ", " << ext_pos.z << ")\n";
    
    for (int iter = 0; iter < max_iter; iter++) {
        // Harmonic force
        Vec3 dr = d_pos - o_pos;
        Vec3 f_harmonic = dr * (-k_drude);
        
        // Coulomb force from external charge
        Vec3 r_ext = ext_pos - d_pos;
        double r_ext_mag = r_ext.norm();
        double f_coulomb_mag = ONE_4PI_EPS0 * qD * ext_charge / (r_ext_mag * r_ext_mag);
        Vec3 f_coulomb = r_ext * (f_coulomb_mag / r_ext_mag);
        
        Vec3 f_total = f_harmonic + f_coulomb;
        double f_mag = f_total.norm();
        
        if (iter % 10 == 0) {
            std::cout << "Iter " << iter << ": |F| = " << f_mag << " kJ/mol/nm\n";
        }
        
        if (f_mag < tolerance) {
            std::cout << "Converged after " << iter << " iterations\n";
            break;
        }
        
        // Update
        double damping = (f_mag > 10.0 * tolerance) ? 0.1 : 0.5;
        Vec3 delta = f_total * (damping / k_drude);
        d_pos = d_pos + delta;
    }
    
    std::cout << "Final Drude displacement from parent: " << (d_pos - o_pos).norm() << " nm\n";
    
    // Calculate induced dipole
    Vec3 dipole = (d_pos - o_pos) * qD;
    std::cout << "Induced dipole magnitude: " << dipole.norm() << " e·nm\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    debugSingleWater();
    return 0;
}