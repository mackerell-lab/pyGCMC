// Correct implementation of SWM4-NDP water with proper Drude SCF theory

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <random>
#include <iomanip>
#include <algorithm>

// Physical constants
const double ONE_4PI_EPS0_KCAL = 332.0637;     // kcal*Å/mol/e²
const double ONE_4PI_EPS0_KJ = 138.935456;     // kJ*nm/mol/e²
const double ANGSTROM_TO_NM = 0.1;
const double KCAL_TO_KJ = 4.184;

// TIP3P parameters for comparison
namespace TIP3P {
    const double qO = -0.834;   
    const double qH = 0.417;    
    const double sigma_O = 3.15061;  // Å
    const double eps_O = 0.6364;     // kJ/mol
    const double rOH = 0.9572;       // Å
    const double aHOH = 104.52;      // degrees
}

// SWM4-NDP parameters (Lamoureux et al. 2003)
namespace SWM4_NDP {
    const double qO = 1.71636;      // Parent oxygen (positive because charge is on Drude)
    const double qD = -1.71636;     // Drude particle
    const double qH = 0.55733;      // Hydrogen
    const double qM = -1.11466;     // Virtual M-site
    const double sigma_O = 3.18395; // Å
    const double eps_O = 0.21094;   // kcal/mol
    const double alpha = 0.978;     // Å³ (polarizability)
    const double k_drude = 1000.0;  // kcal/mol/Å² (force constant)
    const double rOH = 0.9572;      // Å
    const double aHOH = 104.52;     // degrees
    const double thole = 2.6;       // Thole screening parameter (if needed)
    
    // M-site weights
    const double w_O = 0.786646558;
    const double w_H = 0.106676721;
}

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

class WaterBox {
protected:
    std::vector<Vec3> positions;
    std::vector<double> charges;
    std::vector<int> types;  
    std::vector<int> molecule_id;  // Which molecule each atom belongs to
    double box_size;
    double cutoff;
    int n_waters;
    int atoms_per_water;
    
public:
    WaterBox(double box_, double cut) : box_size(box_), cutoff(cut) {}
    virtual ~WaterBox() = default;
    
    Vec3 minimumImage(const Vec3& r) const {
        Vec3 result = r;
        if (result.x > box_size/2) result.x -= box_size;
        else if (result.x < -box_size/2) result.x += box_size;
        if (result.y > box_size/2) result.y -= box_size;
        else if (result.y < -box_size/2) result.y += box_size;
        if (result.z > box_size/2) result.z -= box_size;
        else if (result.z < -box_size/2) result.z += box_size;
        return result;
    }
    
    void applyPBC(Vec3& pos) {
        while (pos.x < 0) pos.x += box_size;
        while (pos.x > box_size) pos.x -= box_size;
        while (pos.y < 0) pos.y += box_size;
        while (pos.y > box_size) pos.y -= box_size;
        while (pos.z < 0) pos.z += box_size;
        while (pos.z > box_size) pos.z -= box_size;
    }
    
    virtual void generateWaterBox(int target_waters) = 0;
    virtual double calculateEnergy() = 0;
    virtual void moveWater(int water_idx, const Vec3& displacement) = 0;
};

class TIP3PWaterBox : public WaterBox {
public:
    TIP3PWaterBox(double box, double cut) : WaterBox(box, cut) {
        atoms_per_water = 3;
    }
    
    void generateWaterBox(int target_waters) override {
        n_waters = target_waters;
        positions.clear();
        charges.clear();
        types.clear();
        molecule_id.clear();
        
        int n_per_dim = std::ceil(std::cbrt(n_waters));
        double spacing = box_size / n_per_dim;
        
        int water_count = 0;
        for (int i = 0; i < n_per_dim && water_count < n_waters; i++) {
            for (int j = 0; j < n_per_dim && water_count < n_waters; j++) {
                for (int k = 0; k < n_per_dim && water_count < n_waters; k++) {
                    Vec3 center((i + 0.5) * spacing,
                               (j + 0.5) * spacing,
                               (k + 0.5) * spacing);
                    
                    addWaterMolecule(center, water_count);
                    water_count++;
                }
            }
        }
    }
    
    void addWaterMolecule(const Vec3& center, int mol_id) {
        // Oxygen
        positions.push_back(center);
        charges.push_back(TIP3P::qO);
        types.push_back(0);
        molecule_id.push_back(mol_id);
        
        // Hydrogen 1
        double angle_rad = TIP3P::aHOH * M_PI / 180.0;
        Vec3 h1 = center;
        h1.x += TIP3P::rOH * ANGSTROM_TO_NM;
        positions.push_back(h1);
        charges.push_back(TIP3P::qH);
        types.push_back(1);
        molecule_id.push_back(mol_id);
        
        // Hydrogen 2
        Vec3 h2 = center;
        h2.x += TIP3P::rOH * ANGSTROM_TO_NM * std::cos(angle_rad);
        h2.y += TIP3P::rOH * ANGSTROM_TO_NM * std::sin(angle_rad);
        positions.push_back(h2);
        charges.push_back(TIP3P::qH);
        types.push_back(1);
        molecule_id.push_back(mol_id);
    }
    
    double calculateEnergy() override {
        double energy = 0.0;
        double cutoff2 = cutoff * cutoff;
        
        for (size_t i = 0; i < positions.size(); i++) {
            for (size_t j = i + 1; j < positions.size(); j++) {
                // Skip intramolecular
                if (molecule_id[i] == molecule_id[j]) continue;
                
                Vec3 rij = minimumImage(positions[j] - positions[i]);
                double r2 = rij.norm2();
                
                if (r2 < cutoff2 && r2 > 1e-10) {
                    double r = std::sqrt(r2);
                    
                    // Coulomb
                    energy += ONE_4PI_EPS0_KJ * charges[i] * charges[j] / r;
                    
                    // LJ (only O-O)
                    if (types[i] == 0 && types[j] == 0) {
                        double sigma_nm = TIP3P::sigma_O * ANGSTROM_TO_NM;
                        double r6 = r2 * r2 * r2;
                        double sigma6 = sigma_nm * sigma_nm * sigma_nm * sigma_nm * sigma_nm * sigma_nm;
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
            positions[start + i] += displacement;
            applyPBC(positions[start + i]);
        }
    }
};

class SWM4NDPWaterBox : public WaterBox {
private:
    std::vector<int> drude_indices;
    std::vector<int> parent_indices;
    
public:
    SWM4NDPWaterBox(double box, double cut) : WaterBox(box, cut) {
        atoms_per_water = 5;  // O, D, H1, H2, M
    }
    
    void generateWaterBox(int target_waters) override {
        n_waters = target_waters;
        positions.clear();
        charges.clear();
        types.clear();
        molecule_id.clear();
        drude_indices.clear();
        parent_indices.clear();
        
        int n_per_dim = std::ceil(std::cbrt(n_waters));
        double spacing = box_size / n_per_dim;
        
        int water_count = 0;
        for (int i = 0; i < n_per_dim && water_count < n_waters; i++) {
            for (int j = 0; j < n_per_dim && water_count < n_waters; j++) {
                for (int k = 0; k < n_per_dim && water_count < n_waters; k++) {
                    Vec3 center((i + 0.5) * spacing,
                               (j + 0.5) * spacing,
                               (k + 0.5) * spacing);
                    
                    addWaterMolecule(center, water_count);
                    water_count++;
                }
            }
        }
        
        // Initialize Drude positions with SCF
        std::cout << "  Initializing Drude positions with SCF...\n";
        runSCF(0.1, 200);  // Tight tolerance for initialization
    }
    
    void addWaterMolecule(const Vec3& center, int mol_id) {
        int start_idx = positions.size();
        
        // Oxygen (parent)
        positions.push_back(center);
        charges.push_back(SWM4_NDP::qO);
        types.push_back(0);
        molecule_id.push_back(mol_id);
        
        // Drude particle - initially at parent position
        positions.push_back(center);
        charges.push_back(SWM4_NDP::qD);
        types.push_back(1);
        molecule_id.push_back(mol_id);
        drude_indices.push_back(start_idx + 1);
        parent_indices.push_back(start_idx);
        
        // Hydrogen 1
        double angle_rad = SWM4_NDP::aHOH * M_PI / 180.0;
        Vec3 h1 = center;
        h1.x += SWM4_NDP::rOH * ANGSTROM_TO_NM;
        positions.push_back(h1);
        charges.push_back(SWM4_NDP::qH);
        types.push_back(2);
        molecule_id.push_back(mol_id);
        
        // Hydrogen 2
        Vec3 h2 = center;
        h2.x += SWM4_NDP::rOH * ANGSTROM_TO_NM * std::cos(angle_rad);
        h2.y += SWM4_NDP::rOH * ANGSTROM_TO_NM * std::sin(angle_rad);
        positions.push_back(h2);
        charges.push_back(SWM4_NDP::qH);
        types.push_back(2);
        molecule_id.push_back(mol_id);
        
        // Virtual M-site
        updateMSite(start_idx);
    }
    
    void updateMSite(int mol_start_idx) {
        // M-site position is weighted average of O, H1, H2
        Vec3 O = positions[mol_start_idx];
        Vec3 H1 = positions[mol_start_idx + 2];
        Vec3 H2 = positions[mol_start_idx + 3];
        
        Vec3 M;
        M.x = SWM4_NDP::w_O * O.x + SWM4_NDP::w_H * H1.x + SWM4_NDP::w_H * H2.x;
        M.y = SWM4_NDP::w_O * O.y + SWM4_NDP::w_H * H1.y + SWM4_NDP::w_H * H2.y;
        M.z = SWM4_NDP::w_O * O.z + SWM4_NDP::w_H * H1.z + SWM4_NDP::w_H * H2.z;
        
        if (mol_start_idx + 4 < positions.size()) {
            positions[mol_start_idx + 4] = M;
        } else {
            positions.push_back(M);
            charges.push_back(SWM4_NDP::qM);
            types.push_back(3);
            molecule_id.push_back(molecule_id[mol_start_idx]);
        }
    }
    
    void runSCF(double tolerance_kjmol_nm = 1.0, int max_iter = 100) {
        // Convert force constant to kJ/mol/nm²
        double k_drude_kj_nm2 = SWM4_NDP::k_drude * KCAL_TO_KJ * 100.0;
        
        std::vector<Vec3> forces(drude_indices.size());
        std::vector<Vec3> old_positions(drude_indices.size());
        
        for (int iter = 0; iter < max_iter; iter++) {
            double max_force = 0.0;
            
            // Save old positions
            for (size_t d = 0; d < drude_indices.size(); d++) {
                old_positions[d] = positions[drude_indices[d]];
            }
            
            // Calculate forces on all Drudes
            for (size_t d = 0; d < drude_indices.size(); d++) {
                int di = drude_indices[d];
                int pi = parent_indices[d];
                
                Vec3 force(0, 0, 0);
                
                // Harmonic restoring force
                Vec3 dr = positions[di] - positions[pi];
                force = dr * (-k_drude_kj_nm2);
                
                // Electrostatic forces from all other atoms
                for (size_t j = 0; j < positions.size(); j++) {
                    if (j == di) continue;  // Skip self
                    
                    // For Drude, we include intramolecular interactions except with parent
                    if (j == pi) continue;  // Skip parent (already in harmonic)
                    
                    Vec3 rij = minimumImage(positions[j] - positions[di]);
                    double r2 = rij.norm2();
                    
                    if (r2 > 1e-10 && r2 < cutoff * cutoff) {
                        double r = std::sqrt(r2);
                        double f_mag = ONE_4PI_EPS0_KJ * charges[di] * charges[j] / (r2 * r);
                        force += rij * (f_mag / r);  // Note: force in direction of rij
                    }
                }
                
                forces[d] = force;
                max_force = std::max(max_force, force.norm());
            }
            
            if (max_force < tolerance_kjmol_nm) {
                if (iter > 0) {
                    // Check convergence by position change
                    double max_displacement = 0.0;
                    for (size_t d = 0; d < drude_indices.size(); d++) {
                        Vec3 delta = positions[drude_indices[d]] - old_positions[d];
                        max_displacement = std::max(max_displacement, delta.norm());
                    }
                    if (max_displacement < 1e-6) {  // 1e-6 nm = 0.00001 Å
                        break;
                    }
                }
            }
            
            // Update positions using force/k as displacement
            for (size_t d = 0; d < drude_indices.size(); d++) {
                double force_mag = forces[d].norm();
                
                if (force_mag > tolerance_kjmol_nm * 0.01) {
                    // Calculate displacement from force
                    Vec3 displacement = forces[d] * (1.0 / k_drude_kj_nm2);
                    
                    // Apply damping for stability
                    double damping = 0.7;  // Standard damping
                    if (force_mag > 1000.0) {
                        damping = 0.1;  // Strong damping for large forces
                    } else if (force_mag > 100.0) {
                        damping = 0.3;
                    }
                    
                    displacement = displacement * damping;
                    
                    // Limit maximum displacement
                    double max_disp = 0.0001;  // 0.001 Å
                    double disp_mag = displacement.norm();
                    if (disp_mag > max_disp) {
                        displacement = displacement * (max_disp / disp_mag);
                    }
                    
                    positions[drude_indices[d]] += displacement;
                }
            }
        }
    }
    
    double calculateEnergy() override {
        // Optimize Drude positions first
        runSCF();
        
        double energy = 0.0;
        double cutoff2 = cutoff * cutoff;
        
        // Convert k_drude to kJ/mol/nm²
        double k_drude_kj_nm2 = SWM4_NDP::k_drude * KCAL_TO_KJ * 100.0;
        
        // 1. Harmonic energy for Drude oscillators
        for (size_t d = 0; d < drude_indices.size(); d++) {
            Vec3 dr = positions[drude_indices[d]] - positions[parent_indices[d]];
            energy += 0.5 * k_drude_kj_nm2 * dr.norm2();
        }
        
        // 2. Pairwise interactions
        for (size_t i = 0; i < positions.size(); i++) {
            for (size_t j = i + 1; j < positions.size(); j++) {
                // Skip intramolecular except for specific cases
                if (molecule_id[i] == molecule_id[j]) {
                    // Check if this is Drude-parent (skip, already in harmonic)
                    bool is_drude_parent = false;
                    for (size_t d = 0; d < drude_indices.size(); d++) {
                        if ((i == drude_indices[d] && j == parent_indices[d]) ||
                            (j == drude_indices[d] && i == parent_indices[d])) {
                            is_drude_parent = true;
                            break;
                        }
                    }
                    if (is_drude_parent) continue;
                    
                    // For SWM4-NDP, we need intramolecular Drude interactions
                    bool has_drude = false;
                    for (size_t d = 0; d < drude_indices.size(); d++) {
                        if (i == drude_indices[d] || j == drude_indices[d]) {
                            has_drude = true;
                            break;
                        }
                    }
                    if (!has_drude) continue;  // Skip other intramolecular
                }
                
                Vec3 rij = minimumImage(positions[j] - positions[i]);
                double r2 = rij.norm2();
                
                if (r2 < cutoff2 && r2 > 1e-10) {
                    double r = std::sqrt(r2);
                    
                    // Coulomb
                    energy += ONE_4PI_EPS0_KJ * charges[i] * charges[j] / r;
                    
                    // LJ (only O-O)
                    if (types[i] == 0 && types[j] == 0) {
                        double sigma_nm = SWM4_NDP::sigma_O * ANGSTROM_TO_NM;
                        double eps_kj = SWM4_NDP::eps_O * KCAL_TO_KJ;
                        double r6 = r2 * r2 * r2;
                        double sigma6 = sigma_nm * sigma_nm * sigma_nm * sigma_nm * sigma_nm * sigma_nm;
                        double r12 = r6 * r6;
                        energy += 4.0 * eps_kj * (sigma6*sigma6/r12 - sigma6/r6);
                    }
                }
            }
        }
        
        return energy;
    }
    
    void moveWater(int water_idx, const Vec3& displacement) override {
        int start = water_idx * atoms_per_water;
        
        // Move O, D, H1, H2 (not M)
        for (int i = 0; i < 4; i++) {
            positions[start + i] += displacement;
            applyPBC(positions[start + i]);
        }
        
        // Update M-site position based on new O, H1, H2 positions
        updateMSite(start);
        
        // M-site might also need PBC
        applyPBC(positions[start + 4]);
    }
};

void runBenchmark() {
    // Parameters
    double box_size_A = 10.0;  // Å
    double box_size_nm = box_size_A * ANGSTROM_TO_NM;
    double cutoff_nm = 0.9;  // 9 Å
    int n_waters = 10;
    int n_moves = 100;
    double max_displacement = 0.02;  // nm
    
    std::cout << "=================================================================\n";
    std::cout << "CHARMM Standard vs CHARMM Drude Performance Benchmark\n";
    std::cout << "=================================================================\n";
    std::cout << "Box size: " << box_size_A << " Å (" << n_waters << " waters)\n";
    std::cout << "Cutoff: " << cutoff_nm/ANGSTROM_TO_NM << " Å\n";
    std::cout << "Number of moves: " << n_moves << "\n\n";
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(42);  // Fixed seed for reproducibility
    std::uniform_real_distribution<> dis(-max_displacement, max_displacement);
    std::uniform_int_distribution<> water_dis(0, n_waters-1);
    
    // Test TIP3P
    std::cout << "Testing TIP3P (Standard CHARMM)...\n";
    TIP3PWaterBox tip3p_box(box_size_nm, cutoff_nm);
    tip3p_box.generateWaterBox(n_waters);
    
    double tip3p_initial = tip3p_box.calculateEnergy();
    std::cout << "  Initial energy: " << tip3p_initial << " kJ/mol\n";
    std::cout << "  Per water: " << tip3p_initial/n_waters << " kJ/mol\n";
    
    auto start = std::chrono::high_resolution_clock::now();
    double tip3p_energy = tip3p_initial;
    
    for (int move = 0; move < n_moves; move++) {
        int water_idx = water_dis(gen);
        Vec3 displacement(dis(gen), dis(gen), dis(gen));
        tip3p_box.moveWater(water_idx, displacement);
        tip3p_energy = tip3p_box.calculateEnergy();
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto tip3p_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    
    std::cout << "  Final energy: " << tip3p_energy << " kJ/mol\n\n";
    
    // Test SWM4-NDP
    std::cout << "Testing SWM4-NDP (CHARMM Drude with SCF)...\n";
    SWM4NDPWaterBox swm4_box(box_size_nm, cutoff_nm);
    swm4_box.generateWaterBox(n_waters);
    
    double swm4_initial = swm4_box.calculateEnergy();
    std::cout << "  Initial energy: " << swm4_initial << " kJ/mol\n";
    std::cout << "  Per water: " << swm4_initial/n_waters << " kJ/mol\n";
    
    // Reset random generator for same moves
    gen.seed(42);
    
    start = std::chrono::high_resolution_clock::now();
    double swm4_energy = swm4_initial;
    
    for (int move = 0; move < n_moves; move++) {
        int water_idx = water_dis(gen);
        Vec3 displacement(dis(gen), dis(gen), dis(gen));
        swm4_box.moveWater(water_idx, displacement);
        swm4_energy = swm4_box.calculateEnergy();
        
        if (move % 20 == 0) {
            std::cout << "  Move " << move << ": E = " << swm4_energy << " kJ/mol\n";
        }
    }
    
    end = std::chrono::high_resolution_clock::now();
    auto swm4_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    
    std::cout << "  Final energy: " << swm4_energy << " kJ/mol\n\n";
    
    // Results
    std::cout << "=================================================================\n";
    std::cout << "RESULTS\n";
    std::cout << "=================================================================\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "TIP3P (Standard CHARMM):\n";
    std::cout << "  Total time: " << tip3p_time/1000.0 << " ms\n";
    std::cout << "  Time per move: " << tip3p_time/1000.0/n_moves << " ms\n";
    std::cout << "  Energy change: " << (tip3p_energy - tip3p_initial) << " kJ/mol\n\n";
    
    std::cout << "SWM4-NDP (CHARMM Drude with SCF):\n";
    std::cout << "  Total time: " << swm4_time/1000.0 << " ms\n";
    std::cout << "  Time per move: " << swm4_time/1000.0/n_moves << " ms\n";
    std::cout << "  Energy change: " << (swm4_energy - swm4_initial) << " kJ/mol\n\n";
    
    std::cout << "Performance ratio (Drude/Standard): " 
              << std::setprecision(1) << (double)swm4_time/tip3p_time << "x slower\n";
}

int main() {
    runBenchmark();
    return 0;
}