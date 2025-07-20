// Debug water box energy calculation

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

void debugTwoWaters() {
    std::cout << "=== Debugging Two SWM4-NDP Water Molecules ===\n\n";
    
    // Parameters
    const double qO = 1.71636;
    const double qD = -1.71636;
    const double qH = 0.55733;
    const double qM = -1.11466;
    const double k_drude = 1000.0 * KCAL_TO_KJ / (ANGSTROM_TO_NM * ANGSTROM_TO_NM);
    const double sigma_O = 0.318395 * ANGSTROM_TO_NM;
    const double eps_O = 0.21094 * KCAL_TO_KJ;
    
    // Create two water molecules 0.5 nm apart
    std::vector<Vec3> pos;
    std::vector<double> charges;
    std::vector<std::string> names;
    std::vector<int> molecule;  // Which molecule each atom belongs to
    std::vector<int> types;     // 0=O, 1=D, 2=H, 3=M
    
    // Water 1 at origin
    pos.push_back(Vec3(0, 0, 0));           charges.push_back(qO);  names.push_back("O1");  molecule.push_back(0); types.push_back(0);
    pos.push_back(Vec3(0, 0, 0.0001));      charges.push_back(qD);  names.push_back("D1");  molecule.push_back(0); types.push_back(1);
    pos.push_back(Vec3(0.09572, 0, 0));     charges.push_back(qH);  names.push_back("H11"); molecule.push_back(0); types.push_back(2);
    pos.push_back(Vec3(-0.024, 0.0927, 0)); charges.push_back(qH);  names.push_back("H12"); molecule.push_back(0); types.push_back(2);
    pos.push_back(Vec3(0.0237, 0.0237, 0)); charges.push_back(qM);  names.push_back("M1");  molecule.push_back(0); types.push_back(3);
    
    // Water 2 at 0.5 nm distance
    double dx = 0.5;
    pos.push_back(Vec3(dx, 0, 0));           charges.push_back(qO);  names.push_back("O2");  molecule.push_back(1); types.push_back(0);
    pos.push_back(Vec3(dx, 0, 0.0001));      charges.push_back(qD);  names.push_back("D2");  molecule.push_back(1); types.push_back(1);
    pos.push_back(Vec3(dx+0.09572, 0, 0));   charges.push_back(qH);  names.push_back("H21"); molecule.push_back(1); types.push_back(2);
    pos.push_back(Vec3(dx-0.024, 0.0927, 0));charges.push_back(qH);  names.push_back("H22"); molecule.push_back(1); types.push_back(2);
    pos.push_back(Vec3(dx+0.0237, 0.0237, 0));charges.push_back(qM);  names.push_back("M2");  molecule.push_back(1); types.push_back(3);
    
    // First, run SCF to optimize Drude positions
    std::cout << "Initial Drude positions:\n";
    std::cout << "  D1: " << pos[1].z << " nm from O1\n";
    std::cout << "  D2: " << pos[6].z << " nm from O2\n\n";
    
    std::cout << "Running SCF...\n";
    
    double tolerance = 0.1;  // kJ/mol/nm
    int max_iter = 100;
    std::vector<int> drude_indices = {1, 6};
    std::vector<int> parent_indices = {0, 5};
    
    for (int iter = 0; iter < max_iter; iter++) {
        double max_force = 0.0;
        
        for (size_t d = 0; d < drude_indices.size(); d++) {
            int di = drude_indices[d];
            int pi = parent_indices[d];
            
            Vec3 force(0, 0, 0);
            
            // Harmonic force
            Vec3 dr = pos[di] - pos[pi];
            force = force - (dr * k_drude);
            
            // Coulomb forces from all other atoms
            for (size_t j = 0; j < pos.size(); j++) {
                if (j == di || j == pi) continue;
                
                // Skip intramolecular for Drude (except we need to include H and M of same molecule!)
                // Actually, in SWM4-NDP, Drude interacts with all atoms except its parent
                
                Vec3 rij = pos[j] - pos[di];
                double r2 = rij.norm2();
                double r = std::sqrt(r2);
                
                if (r > 1e-10) {  // Avoid division by zero
                    double f_mag = ONE_4PI_EPS0 * charges[di] * charges[j] / (r2 * r);
                    force = force - (rij * (f_mag / r));
                }
            }
            
            double f_mag = force.norm();
            max_force = std::max(max_force, f_mag);
            
            // Update
            if (f_mag > tolerance * 0.01) {
                double damping = (f_mag > 10.0 * tolerance) ? 0.1 : 0.5;
                Vec3 delta = force * (damping / k_drude);
                pos[di] = pos[di] + delta;
            }
        }
        
        if (iter % 20 == 0) {
            std::cout << "  Iter " << iter << ": max force = " << max_force << " kJ/mol/nm\n";
        }
        
        if (max_force < tolerance) {
            std::cout << "  Converged after " << iter << " iterations\n";
            break;
        }
    }
    
    std::cout << "\nFinal Drude displacements:\n";
    std::cout << "  D1: " << (pos[1] - pos[0]).norm() << " nm from O1\n";
    std::cout << "  D2: " << (pos[6] - pos[5]).norm() << " nm from O2\n\n";
    
    // Now calculate energy
    std::cout << "Calculating energy components:\n\n";
    
    double e_harmonic = 0.0;
    double e_coulomb = 0.0;
    double e_lj = 0.0;
    
    // Harmonic energy
    for (size_t d = 0; d < drude_indices.size(); d++) {
        Vec3 dr = pos[drude_indices[d]] - pos[parent_indices[d]];
        double e = 0.5 * k_drude * dr.norm2();
        e_harmonic += e;
        std::cout << "  Harmonic " << d+1 << ": " << e << " kJ/mol\n";
    }
    
    std::cout << "\nPairwise interactions:\n";
    
    // Pairwise interactions
    for (size_t i = 0; i < pos.size(); i++) {
        for (size_t j = i + 1; j < pos.size(); j++) {
            // Skip intramolecular
            if (molecule[i] == molecule[j]) {
                // But we need to check if this is a Drude-parent pair
                bool is_drude_parent = false;
                for (size_t d = 0; d < drude_indices.size(); d++) {
                    if ((i == drude_indices[d] && j == parent_indices[d]) ||
                        (j == drude_indices[d] && i == parent_indices[d])) {
                        is_drude_parent = true;
                        break;
                    }
                }
                
                if (!is_drude_parent) {
                    continue;  // Skip other intramolecular
                } else {
                    std::cout << "  Note: Skipping Drude-parent Coulomb for " << names[i] << "-" << names[j] << "\n";
                    continue;  // Skip Drude-parent Coulomb (already in harmonic)
                }
            }
            
            Vec3 rij = pos[j] - pos[i];
            double r = rij.norm();
            
            // Coulomb
            double e_coul = ONE_4PI_EPS0 * charges[i] * charges[j] / r;
            e_coulomb += e_coul;
            
            if (std::abs(e_coul) > 0.1) {  // Only print significant interactions
                std::cout << "  Coulomb " << names[i] << "-" << names[j] 
                          << ": " << e_coul << " kJ/mol (r=" << r << " nm)\n";
            }
            
            // LJ (only O-O)
            if (types[i] == 0 && types[j] == 0) {
                double r2 = r * r;
                double sigma6 = sigma_O * sigma_O * sigma_O * sigma_O * sigma_O * sigma_O;
                double r6 = r2 * r2 * r2;
                double r12 = r6 * r6;
                double e_lj_pair = 4.0 * eps_O * (sigma6*sigma6/r12 - sigma6/r6);
                e_lj += e_lj_pair;
                std::cout << "  LJ " << names[i] << "-" << names[j] 
                          << ": " << e_lj_pair << " kJ/mol\n";
            }
        }
    }
    
    std::cout << "\nEnergy summary:\n";
    std::cout << "  Harmonic: " << e_harmonic << " kJ/mol\n";
    std::cout << "  Coulomb:  " << e_coulomb << " kJ/mol\n";
    std::cout << "  LJ:       " << e_lj << " kJ/mol\n";
    std::cout << "  Total:    " << (e_harmonic + e_coulomb + e_lj) << " kJ/mol\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    debugTwoWaters();
    return 0;
}