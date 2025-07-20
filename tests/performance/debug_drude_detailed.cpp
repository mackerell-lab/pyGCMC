// Detailed debugging of Drude implementation

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

const double ONE_4PI_EPS0_KCAL = 332.0637;
const double ONE_4PI_EPS0_KJ = 138.935456;
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

void debugSingleWaterDetailed() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=== Detailed Debug of Single SWM4-NDP Water ===\n\n";
    
    // Parameters
    const double qO = 1.71636;
    const double qD = -1.71636;
    const double qH = 0.55733;
    const double qM = -1.11466;
    const double k_drude_kcal = 1000.0;  // kcal/mol/Å²
    const double k_drude_kj_nm2 = k_drude_kcal * KCAL_TO_KJ * 100.0;
    
    std::cout << "Parameters:\n";
    std::cout << "  k_drude = " << k_drude_kcal << " kcal/mol/Å²\n";
    std::cout << "  k_drude = " << k_drude_kj_nm2 << " kJ/mol/nm²\n";
    std::cout << "  Total charge = " << (qO + qD + 2*qH + qM) << "\n\n";
    
    // Create single water
    std::vector<Vec3> pos;
    std::vector<double> charges;
    std::vector<std::string> names;
    
    // Positions in nm
    pos.push_back(Vec3(0, 0, 0));           charges.push_back(qO);  names.push_back("O");
    pos.push_back(Vec3(0, 0, 0));           charges.push_back(qD);  names.push_back("D");
    pos.push_back(Vec3(0.09572, 0, 0));     charges.push_back(qH);  names.push_back("H1");
    pos.push_back(Vec3(-0.024, 0.0927, 0)); charges.push_back(qH);  names.push_back("H2");
    
    // Calculate M position
    Vec3 M;
    M.x = 0.786646558 * pos[0].x + 0.106676721 * pos[2].x + 0.106676721 * pos[3].x;
    M.y = 0.786646558 * pos[0].y + 0.106676721 * pos[2].y + 0.106676721 * pos[3].y;
    M.z = 0.786646558 * pos[0].z + 0.106676721 * pos[2].z + 0.106676721 * pos[3].z;
    pos.push_back(M); charges.push_back(qM); names.push_back("M");
    
    std::cout << "Initial positions (nm):\n";
    for (size_t i = 0; i < pos.size(); i++) {
        std::cout << "  " << names[i] << ": (" << pos[i].x << ", " << pos[i].y << ", " << pos[i].z << ")\n";
    }
    
    // Calculate intramolecular energy (no SCF yet)
    std::cout << "\nIntramolecular energies (Drude at parent):\n";
    double e_total = 0.0;
    
    // Harmonic
    Vec3 dr = pos[1] - pos[0];
    double e_harm = 0.5 * k_drude_kj_nm2 * dr.norm2();
    std::cout << "  Harmonic O-D: " << e_harm << " kJ/mol\n";
    e_total += e_harm;
    
    // Coulomb (skip O-D, it's in harmonic)
    for (size_t i = 0; i < pos.size(); i++) {
        for (size_t j = i + 1; j < pos.size(); j++) {
            if (i == 0 && j == 1) continue;  // Skip O-D
            
            Vec3 rij = pos[j] - pos[i];
            double r = rij.norm();
            double e_coul = ONE_4PI_EPS0_KJ * charges[i] * charges[j] / r;
            
            std::cout << "  Coulomb " << names[i] << "-" << names[j] 
                      << ": " << e_coul << " kJ/mol (r=" << r << " nm)\n";
            e_total += e_coul;
        }
    }
    
    std::cout << "  Total intramolecular: " << e_total << " kJ/mol\n";
    
    // Now do SCF
    std::cout << "\nSCF optimization:\n";
    double tolerance = 0.1;  // kJ/mol/nm
    
    for (int iter = 0; iter < 50; iter++) {
        // Calculate force on Drude
        Vec3 force(0, 0, 0);
        
        // Harmonic force
        dr = pos[1] - pos[0];
        Vec3 f_harm = dr * (-k_drude_kj_nm2);
        force = force + f_harm;
        
        // Coulomb forces from H1, H2, M
        for (int j = 2; j < 5; j++) {
            Vec3 rij = pos[j] - pos[1];
            double r = rij.norm();
            double f_mag = ONE_4PI_EPS0_KJ * charges[1] * charges[j] / (r * r);
            force = force + rij * (f_mag / r);
        }
        
        double f_total = force.norm();
        
        if (iter % 10 == 0) {
            std::cout << "  Iter " << iter << ": |F| = " << f_total 
                      << " kJ/mol/nm, r_OD = " << dr.norm() << " nm\n";
        }
        
        if (f_total < tolerance) {
            std::cout << "  Converged at iter " << iter << "\n";
            break;
        }
        
        // Update Drude position
        Vec3 displacement = force * (1.0 / k_drude_kj_nm2);
        double damping = (f_total > 100.0) ? 0.1 : 0.5;
        displacement = displacement * damping;
        
        // Limit displacement
        double max_disp = 0.0001;  // nm
        double disp_mag = displacement.norm();
        if (disp_mag > max_disp) {
            displacement = displacement * (max_disp / disp_mag);
        }
        
        pos[1] = pos[1] + displacement;
    }
    
    // Final energy
    std::cout << "\nFinal energy after SCF:\n";
    e_total = 0.0;
    
    // Harmonic
    dr = pos[1] - pos[0];
    e_harm = 0.5 * k_drude_kj_nm2 * dr.norm2();
    std::cout << "  Harmonic O-D: " << e_harm << " kJ/mol\n";
    std::cout << "  O-D distance: " << dr.norm() << " nm (" << dr.norm()/ANGSTROM_TO_NM << " Å)\n";
    e_total += e_harm;
    
    // All Coulomb
    for (size_t i = 0; i < pos.size(); i++) {
        for (size_t j = i + 1; j < pos.size(); j++) {
            if (i == 0 && j == 1) continue;  // Skip O-D
            
            Vec3 rij = pos[j] - pos[i];
            double r = rij.norm();
            double e_coul = ONE_4PI_EPS0_KJ * charges[i] * charges[j] / r;
            
            if (std::abs(e_coul) > 0.1) {
                std::cout << "  Coulomb " << names[i] << "-" << names[j] 
                          << ": " << e_coul << " kJ/mol\n";
            }
            e_total += e_coul;
        }
    }
    
    std::cout << "\n  Total energy: " << e_total << " kJ/mol\n";
    
    // Induced dipole
    Vec3 dipole = (pos[1] - pos[0]) * qD;
    std::cout << "\nInduced dipole: " << dipole.norm() << " e·nm\n";
}

int main() {
    debugSingleWaterDetailed();
    return 0;
}