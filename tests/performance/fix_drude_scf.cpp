// Fix the SCF algorithm for Drude

#include <iostream>
#include <cmath>
#include <iomanip>

const double ONE_4PI_EPS0_KJ = 138.935456;     // kJ*nm/mol/e²
const double KCAL_TO_KJ = 4.184;

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }
    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator*(double s) const { return Vec3(x * s, y * s, z * s); }
    Vec3& operator+=(const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    double norm2() const { return x*x + y*y + z*z; }
};

void debugForceCalculation() {
    std::cout << "=== Debugging Force Calculation ===\n\n";
    
    // Parameters
    const double qO = 1.71636;
    const double qD = -1.71636;
    const double qH = 0.55733;
    const double qM = -1.11466;
    const double k_drude = 1000.0 * KCAL_TO_KJ * 100.0;  // kJ/mol/nm²
    
    // Simple case: Drude at origin, others at unit positions
    Vec3 D(0, 0, 0);
    Vec3 O(0, 0, 0);
    Vec3 H1(0.1, 0, 0);     // 0.1 nm = 1 Å
    Vec3 H2(-0.025, 0.095, 0);
    Vec3 M(0.008, 0.010, 0);
    
    std::cout << "Test configuration:\n";
    std::cout << "  D at origin\n";
    std::cout << "  O at origin (parent)\n";
    std::cout << "  H1 at (0.1, 0, 0) nm\n";
    std::cout << "  H2 at (-0.025, 0.095, 0) nm\n";
    std::cout << "  M at (0.008, 0.01, 0) nm\n\n";
    
    // Calculate individual forces
    std::cout << "Individual forces on Drude:\n";
    
    // 1. Harmonic force (restoring to parent)
    Vec3 dr = D - O;
    Vec3 f_harm = dr * (-k_drude);
    std::cout << "  Harmonic: F = -k*dr = 0 (D at O)\n";
    
    // 2. Coulomb from H1
    Vec3 r_DH1 = H1 - D;
    double dist_DH1 = r_DH1.norm();
    double f_mag_DH1 = ONE_4PI_EPS0_KJ * qD * qH / (dist_DH1 * dist_DH1);
    Vec3 f_DH1 = r_DH1 * (f_mag_DH1 / dist_DH1);
    
    std::cout << "  From H1:\n";
    std::cout << "    r_DH1 = (" << r_DH1.x << ", " << r_DH1.y << ", " << r_DH1.z << ") nm\n";
    std::cout << "    |r_DH1| = " << dist_DH1 << " nm\n";
    std::cout << "    qD*qH = " << qD << " * " << qH << " = " << qD*qH << "\n";
    std::cout << "    Force magnitude = " << f_mag_DH1 << " kJ/mol/nm\n";
    std::cout << "    Force direction: toward H1 (attractive)\n";
    std::cout << "    Force vector: (" << f_DH1.x << ", " << f_DH1.y << ", " << f_DH1.z << ")\n";
    
    // Check: negative charges attract positive charges
    std::cout << "\n  Force direction check:\n";
    std::cout << "    D has charge " << qD << " (negative)\n";
    std::cout << "    H has charge " << qH << " (positive)\n";
    std::cout << "    They should ATTRACT\n";
    std::cout << "    Force on D should point TOWARD H\n";
    
    // The issue: In Coulomb's law, F = k*q1*q2/r² * r̂
    // For opposite charges, q1*q2 < 0, so F points opposite to r̂
    // This means attractive force!
    
    std::cout << "\n=== CORRECTED Force Calculation ===\n";
    std::cout << "The force on particle i due to particle j is:\n";
    std::cout << "  F_i = k*qi*qj/r² * r̂_ij\n";
    std::cout << "where r̂_ij points FROM i TO j\n\n";
    
    std::cout << "For D-H interaction:\n";
    std::cout << "  F_D = k*qD*qH/r² * r̂_DH\n";
    std::cout << "  Since qD < 0 and qH > 0, qD*qH < 0\n";
    std::cout << "  So F_D points opposite to r̂_DH (toward H)\n";
    
    // Correct SCF iteration
    std::cout << "\n=== Correct SCF Iteration ===\n";
    
    // Reset positions
    Vec3 D_pos = O;  // Start at parent
    
    for (int iter = 0; iter < 20; iter++) {
        Vec3 f_total(0, 0, 0);
        
        // Harmonic
        Vec3 dr = D_pos - O;
        f_total += dr * (-k_drude);
        
        // Coulomb forces - CORRECTED
        // From H1
        Vec3 r_vec = H1 - D_pos;
        double r = r_vec.norm();
        if (r > 1e-10) {
            // Force magnitude (always positive)
            double f_coulomb = ONE_4PI_EPS0_KJ * std::abs(qD * qH) / (r * r);
            // Direction: attractive for opposite charges
            Vec3 f_dir = r_vec * (1.0/r);  // Unit vector from D to H
            if (qD * qH < 0) {  // Opposite charges - attractive
                f_total += f_dir * f_coulomb;
            } else {  // Same charges - repulsive
                f_total += f_dir * (-f_coulomb);
            }
        }
        
        // From H2
        r_vec = H2 - D_pos;
        r = r_vec.norm();
        if (r > 1e-10) {
            double f_coulomb = ONE_4PI_EPS0_KJ * std::abs(qD * qH) / (r * r);
            Vec3 f_dir = r_vec * (1.0/r);
            if (qD * qH < 0) {
                f_total += f_dir * f_coulomb;
            } else {
                f_total += f_dir * (-f_coulomb);
            }
        }
        
        // From M
        r_vec = M - D_pos;
        r = r_vec.norm();
        if (r > 1e-10) {
            double f_coulomb = ONE_4PI_EPS0_KJ * std::abs(qD * qM) / (r * r);
            Vec3 f_dir = r_vec * (1.0/r);
            if (qD * qM < 0) {  // Same sign - repulsive
                f_total += f_dir * (-f_coulomb);
            } else {
                f_total += f_dir * f_coulomb;
            }
        }
        
        double f_mag = f_total.norm();
        
        if (iter % 5 == 0) {
            std::cout << "Iter " << iter << ": |F| = " << f_mag << " kJ/mol/nm, "
                      << "D at (" << D_pos.x*1000 << ", " << D_pos.y*1000 << ", " 
                      << D_pos.z*1000 << ") pm\n";
        }
        
        if (f_mag < 1.0) break;
        
        // Update: x_new = x_old + F/k
        Vec3 displacement = f_total * (1.0 / k_drude);
        D_pos += displacement * 0.5;  // Damping
    }
    
    std::cout << "\nFinal Drude position: (" << D_pos.x*1000 << ", " << D_pos.y*1000 
              << ", " << D_pos.z*1000 << ") pm from O\n";
    std::cout << "Displacement: " << (D_pos - O).norm()*1000 << " pm\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    debugForceCalculation();
    return 0;
}