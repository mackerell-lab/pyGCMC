// Deep analysis of Drude physics in SWM4-NDP

#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

const double ONE_4PI_EPS0_KJ = 138.935456;     // kJ*nm/mol/e²
const double ANGSTROM_TO_NM = 0.1;
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
    void normalize() { double n = norm(); if (n > 0) { x /= n; y /= n; z /= n; } }
};

void analyzeIdealDrudePosition() {
    std::cout << "=== Analyzing Ideal Drude Position in SWM4-NDP ===\n\n";
    
    // Parameters
    const double qO = 1.71636;
    const double qD = -1.71636;
    const double qH = 0.55733;
    const double qM = -1.11466;
    const double k_drude = 1000.0;  // kcal/mol/Å²
    const double k_drude_kj_nm2 = k_drude * KCAL_TO_KJ * 100.0;
    
    // Water geometry (in nm)
    const double rOH = 0.09572;  // nm
    const double aHOH = 104.52 * M_PI / 180.0;
    
    // Build water molecule at origin
    Vec3 O(0, 0, 0);
    Vec3 H1(rOH, 0, 0);
    Vec3 H2(rOH * std::cos(aHOH), rOH * std::sin(aHOH), 0);
    
    // Calculate M position
    const double w_O = 0.786646558;
    const double w_H = 0.106676721;
    Vec3 M = O * w_O + H1 * w_H + H2 * w_H;
    
    std::cout << "Water geometry (nm):\n";
    std::cout << "  O:  (" << O.x << ", " << O.y << ", " << O.z << ")\n";
    std::cout << "  H1: (" << H1.x << ", " << H1.y << ", " << H1.z << ")\n";
    std::cout << "  H2: (" << H2.x << ", " << H2.y << ", " << H2.z << ")\n";
    std::cout << "  M:  (" << M.x << ", " << M.y << ", " << M.z << ")\n\n";
    
    // Analyze forces at different Drude positions
    std::cout << "Force analysis at different Drude positions:\n\n";
    
    // Test positions along different directions
    std::vector<Vec3> test_directions = {
        Vec3(1, 0, 0),      // Along x
        Vec3(0, 1, 0),      // Along y
        Vec3(0, 0, 1),      // Along z
        Vec3(1, 1, 0),      // xy diagonal
        (H1 + H2) * 0.5 - O // Toward H atoms
    };
    
    for (auto& dir : test_directions) {
        dir.normalize();
    }
    
    std::vector<double> displacements = {0.0, 0.0001, 0.0005, 0.001, 0.002, 0.005};  // nm
    
    for (size_t d = 0; d < test_directions.size(); d++) {
        std::cout << "Direction " << d+1 << ": (" << test_directions[d].x 
                  << ", " << test_directions[d].y << ", " << test_directions[d].z << ")\n";
        
        for (double disp : displacements) {
            Vec3 D = O + test_directions[d] * disp;
            
            // Calculate forces
            Vec3 f_total(0, 0, 0);
            
            // 1. Harmonic force (toward O)
            Vec3 dr = D - O;
            Vec3 f_harm = dr * (-k_drude_kj_nm2);
            f_total += f_harm;
            
            // 2. Coulomb from H1
            Vec3 r_DH1 = H1 - D;
            double dist_DH1 = r_DH1.norm();
            if (dist_DH1 > 1e-10) {
                double f_mag = ONE_4PI_EPS0_KJ * qD * qH / (dist_DH1 * dist_DH1);
                Vec3 f_DH1 = r_DH1 * (f_mag / dist_DH1);
                f_total += f_DH1;
            }
            
            // 3. Coulomb from H2
            Vec3 r_DH2 = H2 - D;
            double dist_DH2 = r_DH2.norm();
            if (dist_DH2 > 1e-10) {
                double f_mag = ONE_4PI_EPS0_KJ * qD * qH / (dist_DH2 * dist_DH2);
                Vec3 f_DH2 = r_DH2 * (f_mag / dist_DH2);
                f_total += f_DH2;
            }
            
            // 4. Coulomb from M
            Vec3 r_DM = M - D;
            double dist_DM = r_DM.norm();
            if (dist_DM > 1e-10) {
                double f_mag = ONE_4PI_EPS0_KJ * qD * qM / (dist_DM * dist_DM);
                Vec3 f_DM = r_DM * (f_mag / dist_DM);
                f_total += f_DM;
            }
            
            if (disp == 0.0 || disp == 0.001) {
                std::cout << "  disp=" << disp*1000 << " pm: |F|=" << f_total.norm() 
                          << " kJ/mol/nm (toward " << test_directions[d].x << "," 
                          << test_directions[d].y << "," << test_directions[d].z << ")\n";
            }
        }
        std::cout << "\n";
    }
    
    // Now find equilibrium by simple iteration
    std::cout << "Finding equilibrium position by iteration:\n";
    
    Vec3 D = O;  // Start at parent
    double tolerance = 0.1;  // kJ/mol/nm
    
    for (int iter = 0; iter < 100; iter++) {
        Vec3 f_total(0, 0, 0);
        
        // Harmonic
        Vec3 dr = D - O;
        f_total += dr * (-k_drude_kj_nm2);
        
        // Coulomb from H1
        Vec3 r_DH1 = H1 - D;
        double dist_DH1 = r_DH1.norm();
        if (dist_DH1 > 1e-10) {
            double f_mag = ONE_4PI_EPS0_KJ * qD * qH / (dist_DH1 * dist_DH1);
            f_total += r_DH1 * (f_mag / dist_DH1);
        }
        
        // Coulomb from H2
        Vec3 r_DH2 = H2 - D;
        double dist_DH2 = r_DH2.norm();
        if (dist_DH2 > 1e-10) {
            double f_mag = ONE_4PI_EPS0_KJ * qD * qH / (dist_DH2 * dist_DH2);
            f_total += r_DH2 * (f_mag / dist_DH2);
        }
        
        // Coulomb from M
        Vec3 r_DM = M - D;
        double dist_DM = r_DM.norm();
        if (dist_DM > 1e-10) {
            double f_mag = ONE_4PI_EPS0_KJ * qD * qM / (dist_DM * dist_DM);
            f_total += r_DM * (f_mag / dist_DM);
        }
        
        double f_mag = f_total.norm();
        
        if (iter % 10 == 0) {
            std::cout << "  Iter " << iter << ": |F|=" << f_mag << " kJ/mol/nm, "
                      << "D at (" << D.x << ", " << D.y << ", " << D.z << "), "
                      << "disp=" << (D-O).norm()*1000 << " pm\n";
        }
        
        if (f_mag < tolerance) {
            std::cout << "Converged at iter " << iter << "\n";
            break;
        }
        
        // Update position - direct force/k displacement
        Vec3 displacement = f_total * (1.0 / k_drude_kj_nm2);
        
        // Apply damping for stability
        double damping = 0.3;
        if (f_mag > 1000) damping = 0.01;
        else if (f_mag > 100) damping = 0.1;
        
        D += displacement * damping;
    }
    
    std::cout << "\nFinal Drude position: (" << D.x << ", " << D.y << ", " << D.z << ")\n";
    std::cout << "Displacement from O: " << (D-O).norm()*1000 << " pm = " 
              << (D-O).norm()/ANGSTROM_TO_NM << " Å\n";
    
    // Calculate energy at equilibrium
    double e_total = 0.0;
    
    // Harmonic
    Vec3 dr = D - O;
    double e_harm = 0.5 * k_drude_kj_nm2 * dr.norm2();
    e_total += e_harm;
    
    // D-H1
    double r_DH1 = (H1 - D).norm();
    double e_DH1 = ONE_4PI_EPS0_KJ * qD * qH / r_DH1;
    e_total += e_DH1;
    
    // D-H2
    double r_DH2 = (H2 - D).norm();
    double e_DH2 = ONE_4PI_EPS0_KJ * qD * qH / r_DH2;
    e_total += e_DH2;
    
    // D-M
    double r_DM = (M - D).norm();
    double e_DM = ONE_4PI_EPS0_KJ * qD * qM / r_DM;
    e_total += e_DM;
    
    std::cout << "\nEnergy at equilibrium:\n";
    std::cout << "  Harmonic: " << e_harm << " kJ/mol\n";
    std::cout << "  D-H1: " << e_DH1 << " kJ/mol\n";
    std::cout << "  D-H2: " << e_DH2 << " kJ/mol\n";
    std::cout << "  D-M: " << e_DM << " kJ/mol\n";
    std::cout << "  Total: " << e_total << " kJ/mol\n";
    
    // Analyze what went wrong in the original code
    std::cout << "\n=== Analysis of Original Problem ===\n";
    std::cout << "In the original output, D was at (49.9878, 49.9842, 50.0) relative to O at (50, 50, 50)\n";
    std::cout << "This is a displacement of ~0.02 nm = 0.2 Å, which is too large.\n";
    std::cout << "The equilibrium should be much closer to O (< 0.01 Å).\n";
    std::cout << "\nPossible issues:\n";
    std::cout << "1. Initial SCF might have wrong sign in force calculation\n";
    std::cout << "2. Or the damping/step size is too large\n";
    std::cout << "3. Or there's a unit conversion error\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    analyzeIdealDrudePosition();
    return 0;
}