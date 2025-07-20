// tests/performance/benchmark_drude_vs_standard.cpp
// Benchmark comparing CHARMM standard vs CHARMM Drude potential calculations

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <random>
#include <iomanip>

// Physical constants
const double ONE_4PI_EPS0 = 138.935456;  // kJ/mol·nm·e^-2
const double ANGSTROM_TO_NM = 0.1;
const double KCAL_TO_KJ = 4.184;

// Water model parameters
namespace TIP3P {
    const double qO = -0.834;   // Oxygen charge
    const double qH = 0.417;    // Hydrogen charge
    const double sigma_O = 0.315061 * ANGSTROM_TO_NM;  // nm
    const double eps_O = 0.636386;  // kJ/mol
    const double rOH = 0.09572;  // nm
    const double aHOH = 104.52 * M_PI / 180.0;  // radians
}

namespace SWM4_NDP {
    const double qO = 1.71636;    // Oxygen charge (positive because charge is on Drude)
    const double qD = -1.71636;   // Drude charge
    const double qH = 0.55733;    // Hydrogen charge  
    const double qM = -1.11466;   // Virtual site charge
    const double sigma_O = 0.318395 * ANGSTROM_TO_NM;  // nm
    const double eps_O = 0.21094 * KCAL_TO_KJ;  // kJ/mol
    // Force constant: 1000 kcal/mol/Å^2 = 1000 * 4.184 / 0.01 kJ/mol/nm^2
    const double k_drude = 1000.0 * KCAL_TO_KJ / (ANGSTROM_TO_NM * ANGSTROM_TO_NM);  // kJ/mol/nm^2
    const double polarizability = qD * qD / (ONE_4PI_EPS0 * k_drude);
    const double rOH = 0.09572;  // nm
    const double aHOH = 104.52 * M_PI / 180.0;  // radians
}

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

class WaterBox {
protected:
    std::vector<Vec3> positions;
    std::vector<double> charges;
    std::vector<double> masses;
    std::vector<int> types;  // 0=O, 1=H, 2=M (virtual), 3=D (Drude)
    Vec3 box;
    double cutoff;
    int n_waters;
    int atoms_per_water;
    
public:
    WaterBox(double box_size, double cutoff_) : 
        box(box_size, box_size, box_size), cutoff(cutoff_) {}
    
    virtual ~WaterBox() = default;
    
    void generateWaterBox(int target_waters) {
        n_waters = target_waters;
        
        // Calculate grid spacing
        int n_per_dim = std::ceil(std::cbrt(n_waters));
        double spacing = box.x / n_per_dim;
        
        // Generate water molecules
        int water_count = 0;
        for (int i = 0; i < n_per_dim && water_count < n_waters; i++) {
            for (int j = 0; j < n_per_dim && water_count < n_waters; j++) {
                for (int k = 0; k < n_per_dim && water_count < n_waters; k++) {
                    Vec3 center(i * spacing + spacing/2,
                               j * spacing + spacing/2,
                               k * spacing + spacing/2);
                    
                    addWaterMolecule(center);
                    water_count++;
                }
            }
        }
    }
    
    virtual void addWaterMolecule(const Vec3& center) = 0;
    virtual double calculateEnergy() = 0;
    virtual void moveWater(int water_idx, const Vec3& displacement) = 0;
    
    Vec3 minimumImage(const Vec3& r) const {
        Vec3 result = r;
        if (result.x > box.x/2) result.x -= box.x;
        else if (result.x < -box.x/2) result.x += box.x;
        if (result.y > box.y/2) result.y -= box.y;
        else if (result.y < -box.y/2) result.y += box.y;
        if (result.z > box.z/2) result.z -= box.z;
        else if (result.z < -box.z/2) result.z += box.z;
        return result;
    }
};

class TIP3PWaterBox : public WaterBox {
public:
    TIP3PWaterBox(double box_size, double cutoff) : WaterBox(box_size, cutoff) {
        atoms_per_water = 3;
    }
    
    void addWaterMolecule(const Vec3& center) override {
        // Oxygen
        positions.push_back(center);
        charges.push_back(TIP3P::qO);
        masses.push_back(15.9994);
        types.push_back(0);
        
        // Hydrogen 1
        Vec3 h1_pos = center;
        h1_pos.x += TIP3P::rOH;
        positions.push_back(h1_pos);
        charges.push_back(TIP3P::qH);
        masses.push_back(1.008);
        types.push_back(1);
        
        // Hydrogen 2
        Vec3 h2_pos = center;
        h2_pos.x += TIP3P::rOH * std::cos(TIP3P::aHOH);
        h2_pos.y += TIP3P::rOH * std::sin(TIP3P::aHOH);
        positions.push_back(h2_pos);
        charges.push_back(TIP3P::qH);
        masses.push_back(1.008);
        types.push_back(1);
    }
    
    double calculateEnergy() override {
        double energy = 0.0;
        double cutoff2 = cutoff * cutoff;
        
        // Loop over all atom pairs
        for (size_t i = 0; i < positions.size(); i++) {
            for (size_t j = i + 1; j < positions.size(); j++) {
                // Skip intramolecular interactions
                if (i/3 == j/3) continue;
                
                Vec3 rij = minimumImage(positions[j] - positions[i]);
                double r2 = rij.norm2();
                
                if (r2 < cutoff2) {
                    double r = std::sqrt(r2);
                    
                    // Coulomb
                    energy += ONE_4PI_EPS0 * charges[i] * charges[j] / r;
                    
                    // LJ (only O-O)
                    if (types[i] == 0 && types[j] == 0) {
                        double sigma6 = std::pow(TIP3P::sigma_O, 6);
                        double r6 = r2 * r2 * r2;
                        double r12 = r6 * r6;
                        energy += 4.0 * TIP3P::eps_O * (sigma6*sigma6/r12 - sigma6/r6);
                    }
                }
            }
        }
        
        return energy;
    }
    
    void moveWater(int water_idx, const Vec3& displacement) override {
        int start = water_idx * atoms_per_water;
        for (int i = 0; i < atoms_per_water; i++) {
            positions[start + i] = positions[start + i] + displacement;
            // Apply PBC
            if (positions[start + i].x < 0) positions[start + i].x += box.x;
            if (positions[start + i].x > box.x) positions[start + i].x -= box.x;
            if (positions[start + i].y < 0) positions[start + i].y += box.y;
            if (positions[start + i].y > box.y) positions[start + i].y -= box.y;
            if (positions[start + i].z < 0) positions[start + i].z += box.z;
            if (positions[start + i].z > box.z) positions[start + i].z -= box.z;
        }
    }
};

class SWM4NDPWaterBox : public WaterBox {
private:
    std::vector<int> drude_indices;  // Drude particle indices
    std::vector<int> parent_indices; // Parent atom indices
    
public:
    SWM4NDPWaterBox(double box_size, double cutoff) : WaterBox(box_size, cutoff) {
        atoms_per_water = 5;  // O, D, H1, H2, M
    }
    
    void addWaterMolecule(const Vec3& center) override {
        int start_idx = positions.size();
        
        // Oxygen
        positions.push_back(center);
        charges.push_back(SWM4_NDP::qO);
        masses.push_back(15.6);
        types.push_back(0);
        
        // Drude on oxygen - very small initial displacement
        Vec3 drude_pos = center;
        drude_pos.z += 1e-6;  // Much smaller displacement (1e-6 nm = 0.00001 Å)
        positions.push_back(drude_pos);
        charges.push_back(SWM4_NDP::qD);
        masses.push_back(0.4);
        types.push_back(3);
        drude_indices.push_back(start_idx + 1);
        parent_indices.push_back(start_idx);
        
        // Hydrogen 1
        Vec3 h1_pos = center;
        h1_pos.x += SWM4_NDP::rOH;
        positions.push_back(h1_pos);
        charges.push_back(SWM4_NDP::qH);
        masses.push_back(1.008);
        types.push_back(1);
        
        // Hydrogen 2  
        Vec3 h2_pos = center;
        h2_pos.x += SWM4_NDP::rOH * std::cos(SWM4_NDP::aHOH);
        h2_pos.y += SWM4_NDP::rOH * std::sin(SWM4_NDP::aHOH);
        positions.push_back(h2_pos);
        charges.push_back(SWM4_NDP::qH);
        masses.push_back(1.008);
        types.push_back(1);
        
        // Virtual site M - weighted average position
        Vec3 o_pos = center;
        double w_o = 0.786646558;
        double w_h = 0.106676721;
        Vec3 m_pos;
        m_pos.x = w_o * o_pos.x + w_h * h1_pos.x + w_h * h2_pos.x;
        m_pos.y = w_o * o_pos.y + w_h * h1_pos.y + w_h * h2_pos.y;
        m_pos.z = w_o * o_pos.z + w_h * h1_pos.z + w_h * h2_pos.z;
        positions.push_back(m_pos);
        charges.push_back(SWM4_NDP::qM);
        masses.push_back(0.0);
        types.push_back(2);
    }
    
    void runSCF(double tolerance = 1.0, int max_iter = 50) {  // tolerance in kJ/mol/nm
        // SCF optimization for Drude particles
        for (int iter = 0; iter < max_iter; iter++) {
            double max_force = 0.0;
            
            // Calculate forces on each Drude particle
            for (size_t d = 0; d < drude_indices.size(); d++) {
                int drude_idx = drude_indices[d];
                int parent_idx = parent_indices[d];
                
                Vec3 force(0, 0, 0);
                
                // Harmonic force from parent
                Vec3 dr = positions[drude_idx] - positions[parent_idx];
                force = force - (dr * SWM4_NDP::k_drude);
                
                // Coulomb forces from all other atoms
                for (size_t j = 0; j < positions.size(); j++) {
                    if (j == drude_idx || j == parent_idx) continue;
                    
                    Vec3 rij = minimumImage(positions[j] - positions[drude_idx]);
                    double r2 = rij.norm2();
                    
                    if (r2 < cutoff * cutoff) {
                        double r = std::sqrt(r2);
                        double f_mag = ONE_4PI_EPS0 * charges[drude_idx] * charges[j] / (r2 * r);
                        force = force - (rij * (f_mag / r));
                    }
                }
                
                // Update Drude position
                double force_mag = force.norm();
                max_force = std::max(max_force, force_mag);
                
                // Damped update with displacement limit
                if (force_mag > tolerance * 0.01) {
                    // Adaptive damping based on force magnitude
                    double damping = 0.5;
                    if (force_mag > 100.0 * tolerance) {
                        damping = 0.01;  // Very small damping for large forces
                    } else if (force_mag > 10.0 * tolerance) {
                        damping = 0.1;
                    }
                    
                    Vec3 delta = force * (damping / SWM4_NDP::k_drude);
                    
                    // Limit maximum displacement per iteration
                    double max_disp = 0.001;  // 0.001 nm = 0.01 Å
                    double delta_mag = delta.norm();
                    if (delta_mag > max_disp) {
                        delta = delta * (max_disp / delta_mag);
                    }
                    
                    positions[drude_idx] = positions[drude_idx] + delta;
                }
            }
            
            if (max_force < tolerance) break;
        }
    }
    
    double calculateEnergy() override {
        // First optimize Drude positions
        runSCF();
        
        double energy = 0.0;
        double cutoff2 = cutoff * cutoff;
        
        // Harmonic energy for Drude oscillators
        for (size_t d = 0; d < drude_indices.size(); d++) {
            Vec3 dr = positions[drude_indices[d]] - positions[parent_indices[d]];
            energy += 0.5 * SWM4_NDP::k_drude * dr.norm2();
        }
        
        // Loop over all atom pairs
        for (size_t i = 0; i < positions.size(); i++) {
            for (size_t j = i + 1; j < positions.size(); j++) {
                // Skip intramolecular interactions
                if (i/5 == j/5) continue;
                
                // Skip Drude-parent pairs (already in harmonic)
                bool is_drude_parent = false;
                for (size_t d = 0; d < drude_indices.size(); d++) {
                    if ((i == drude_indices[d] && j == parent_indices[d]) ||
                        (j == drude_indices[d] && i == parent_indices[d])) {
                        is_drude_parent = true;
                        break;
                    }
                }
                if (is_drude_parent) continue;
                
                Vec3 rij = minimumImage(positions[j] - positions[i]);
                double r2 = rij.norm2();
                
                if (r2 < cutoff2) {
                    double r = std::sqrt(r2);
                    
                    // Coulomb
                    energy += ONE_4PI_EPS0 * charges[i] * charges[j] / r;
                    
                    // LJ (only O-O)
                    if (types[i] == 0 && types[j] == 0) {
                        double sigma6 = std::pow(SWM4_NDP::sigma_O, 6);
                        double r6 = r2 * r2 * r2;
                        double r12 = r6 * r6;
                        energy += 4.0 * SWM4_NDP::eps_O * (sigma6*sigma6/r12 - sigma6/r6);
                    }
                }
            }
        }
        
        return energy;
    }
    
    void moveWater(int water_idx, const Vec3& displacement) override {
        int start = water_idx * atoms_per_water;
        
        // Move all atoms of the water molecule
        for (int i = 0; i < atoms_per_water; i++) {
            positions[start + i] = positions[start + i] + displacement;
            // Apply PBC
            if (positions[start + i].x < 0) positions[start + i].x += box.x;
            if (positions[start + i].x > box.x) positions[start + i].x -= box.x;
            if (positions[start + i].y < 0) positions[start + i].y += box.y;
            if (positions[start + i].y > box.y) positions[start + i].y -= box.y;
            if (positions[start + i].z < 0) positions[start + i].z += box.z;
            if (positions[start + i].z > box.z) positions[start + i].z -= box.z;
        }
        
        // After moving, need to re-optimize Drude positions
        // (This is done automatically in calculateEnergy)
    }
};

void runBenchmark() {
    // Parameters
    double box_size = 10.0 * ANGSTROM_TO_NM;  // 100 Å³ = 10×10×10 Å
    double cutoff = 0.9;  // 9 Å cutoff
    int n_waters = 33;  // ~33 waters in 100 Å³ for typical density
    int n_moves = 1000;  // Number of MC moves
    double max_displacement = 0.02;  // Maximum displacement per move (nm)
    
    // Start with fewer waters to test stability
    n_waters = 10;
    n_moves = 100;
    
    std::cout << "=================================================================\n";
    std::cout << "CHARMM Standard vs CHARMM Drude Performance Benchmark\n";
    std::cout << "=================================================================\n";
    std::cout << "Box size: " << box_size/ANGSTROM_TO_NM << " Å (" << n_waters << " waters)\n";
    std::cout << "Cutoff: " << cutoff/ANGSTROM_TO_NM << " Å\n";
    std::cout << "Number of moves: " << n_moves << "\n\n";
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-max_displacement, max_displacement);
    std::uniform_int_distribution<> water_dis(0, n_waters-1);
    
    // Test TIP3P (standard CHARMM)
    std::cout << "Testing TIP3P (Standard CHARMM)...\n";
    TIP3PWaterBox tip3p_box(box_size, cutoff);
    tip3p_box.generateWaterBox(n_waters);
    
    // Initial energy
    double initial_tip3p = tip3p_box.calculateEnergy();
    std::cout << "  Initial energy: " << initial_tip3p << " kJ/mol\n";
    
    auto start = std::chrono::high_resolution_clock::now();
    double tip3p_energy = 0.0;
    
    for (int move = 0; move < n_moves; move++) {
        int water_idx = water_dis(gen);
        Vec3 displacement(dis(gen), dis(gen), dis(gen));
        tip3p_box.moveWater(water_idx, displacement);
        tip3p_energy = tip3p_box.calculateEnergy();
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto tip3p_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    
    // Test SWM4-NDP (CHARMM Drude)
    std::cout << "\nTesting SWM4-NDP (CHARMM Drude with SCF)...\n";
    SWM4NDPWaterBox swm4_box(box_size, cutoff);
    swm4_box.generateWaterBox(n_waters);
    
    // Initial energy
    double initial_swm4 = swm4_box.calculateEnergy();
    std::cout << "  Initial energy: " << initial_swm4 << " kJ/mol\n";
    
    start = std::chrono::high_resolution_clock::now();
    double swm4_energy = 0.0;
    
    for (int move = 0; move < n_moves; move++) {
        int water_idx = water_dis(gen);
        Vec3 displacement(dis(gen), dis(gen), dis(gen));
        swm4_box.moveWater(water_idx, displacement);
        swm4_energy = swm4_box.calculateEnergy();  // Includes SCF
    }
    
    end = std::chrono::high_resolution_clock::now();
    auto swm4_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    
    // Results
    std::cout << "\n=================================================================\n";
    std::cout << "RESULTS\n";
    std::cout << "=================================================================\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "TIP3P (Standard CHARMM):\n";
    std::cout << "  Total time: " << tip3p_time/1000.0 << " ms\n";
    std::cout << "  Time per move: " << tip3p_time/1000.0/n_moves << " ms\n";
    std::cout << "  Final energy: " << tip3p_energy << " kJ/mol\n\n";
    
    std::cout << "SWM4-NDP (CHARMM Drude with SCF):\n";
    std::cout << "  Total time: " << swm4_time/1000.0 << " ms\n";
    std::cout << "  Time per move: " << swm4_time/1000.0/n_moves << " ms\n";
    std::cout << "  Final energy: " << swm4_energy << " kJ/mol\n\n";
    
    std::cout << "Performance ratio (Drude/Standard): " 
              << std::setprecision(1) << (double)swm4_time/tip3p_time << "x slower\n";
    
    std::cout << "\n=================================================================\n";
    std::cout << "ANALYSIS\n";
    std::cout << "=================================================================\n";
    std::cout << "The Drude model is slower due to:\n";
    std::cout << "1. More atoms per water (5 vs 3)\n";
    std::cout << "2. SCF optimization required after each move\n";
    std::cout << "3. Additional harmonic terms for Drude oscillators\n";
    std::cout << "4. More complex force calculations\n";
}

int main() {
    runBenchmark();
    return 0;
}