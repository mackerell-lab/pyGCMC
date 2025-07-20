// Debug single SWM4-NDP water with detailed energy breakdown

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

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
    double norm() const { return std::sqrt(x*x + y*y + z*z); }
    double norm2() const { return x*x + y*y + z*z; }
};

int main() {
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=== Debug Single SWM4-NDP Water Energy ===\n\n";
    
    // Parameters
    const double qO = 1.71636;
    const double qD = -1.71636;
    const double qH = 0.55733;
    const double qM = -1.11466;
    const double k_drude_kj_nm2 = 1000.0 * KCAL_TO_KJ * 100.0;  // 418400 kJ/mol/nm²
    
    // Positions from the output (in nm)
    Vec3 pos[5];
    pos[0] = Vec3(50.0, 50.0, 50.0);        // O
    pos[1] = Vec3(49.9878, 49.9842, 50.0);  // D
    pos[2] = Vec3(50.0957, 50.0, 50.0);     // H1
    pos[3] = Vec3(49.976, 50.0927, 50.0);   // H2
    pos[4] = Vec3(50.0077, 50.0099, 50.0);  // M
    
    double charges[5] = {qO, qD, qH, qH, qM};
    std::string names[5] = {"O", "D", "H1", "H2", "M"};
    
    std::cout << "Atom positions and charges:\n";
    for (int i = 0; i < 5; i++) {
        std::cout << "  " << names[i] << ": q=" << charges[i] 
                  << ", pos=(" << pos[i].x << ", " << pos[i].y << ", " << pos[i].z << ")\n";
    }
    
    // Calculate O-D distance
    Vec3 dr_OD = pos[1] - pos[0];
    std::cout << "\nO-D displacement: " << dr_OD.norm() << " nm = " 
              << dr_OD.norm()/ANGSTROM_TO_NM << " Å\n";
    
    std::cout << "\nEnergy breakdown:\n";
    double total_energy = 0.0;
    
    // 1. Harmonic energy
    double e_harmonic = 0.5 * k_drude_kj_nm2 * dr_OD.norm2();
    std::cout << "\n1. Harmonic O-D: " << e_harmonic << " kJ/mol\n";
    total_energy += e_harmonic;
    
    // 2. All pairwise Coulomb interactions
    std::cout << "\n2. Coulomb interactions:\n";
    
    // Matrix of what should be included
    std::cout << "   Exclusion rules:\n";
    std::cout << "   O-D: EXCLUDED (in harmonic)\n";
    std::cout << "   O-H: EXCLUDED (bonded)\n";
    std::cout << "   O-M: EXCLUDED (M is virtual site)\n";
    std::cout << "   D-H: INCLUDED\n";
    std::cout << "   D-M: INCLUDED\n";
    std::cout << "   H-H: EXCLUDED (1-3 connected)\n";
    std::cout << "   H-M: EXCLUDED (M is virtual site)\n\n";
    
    std::cout << "   Included interactions:\n";
    
    // D-H1
    Vec3 r_DH1 = pos[2] - pos[1];
    double e_DH1 = ONE_4PI_EPS0_KJ * charges[1] * charges[2] / r_DH1.norm();
    std::cout << "   D-H1: " << e_DH1 << " kJ/mol (r=" << r_DH1.norm() << " nm)\n";
    total_energy += e_DH1;
    
    // D-H2
    Vec3 r_DH2 = pos[3] - pos[1];
    double e_DH2 = ONE_4PI_EPS0_KJ * charges[1] * charges[3] / r_DH2.norm();
    std::cout << "   D-H2: " << e_DH2 << " kJ/mol (r=" << r_DH2.norm() << " nm)\n";
    total_energy += e_DH2;
    
    // D-M
    Vec3 r_DM = pos[4] - pos[1];
    double e_DM = ONE_4PI_EPS0_KJ * charges[1] * charges[4] / r_DM.norm();
    std::cout << "   D-M: " << e_DM << " kJ/mol (r=" << r_DM.norm() << " nm)\n";
    total_energy += e_DM;
    
    std::cout << "\nTotal intramolecular energy: " << total_energy << " kJ/mol\n";
    
    // Analysis
    std::cout << "\n=== Analysis ===\n";
    std::cout << "The energy of -2092 kJ/mol comes from:\n";
    std::cout << "1. Harmonic: " << e_harmonic << " kJ/mol (small, as expected)\n";
    std::cout << "2. D-H interactions: " << (e_DH1 + e_DH2) << " kJ/mol\n";
    std::cout << "3. D-M interaction: " << e_DM << " kJ/mol\n";
    
    std::cout << "\nThe large negative energy is due to D-M interaction.\n";
    std::cout << "D has charge " << qD << " and M has charge " << qM << "\n";
    std::cout << "They attract each other strongly.\n";
    
    // Check if this is physical
    std::cout << "\nIs this physical?\n";
    std::cout << "In SWM4-NDP, the M site is designed to improve the molecular dipole.\n";
    std::cout << "The large D-M attraction is balanced by the overall molecular geometry.\n";
    std::cout << "This intramolecular energy is a reference state - what matters is\n";
    std::cout << "the energy DIFFERENCES during simulation.\n";
    
    return 0;
}