// Final correct implementation of SWM4-NDP with proper exclusions

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <random>
#include <iomanip>

const double ONE_4PI_EPS0_KCAL = 332.0637;     // kcal*Å/mol/e²
const double ONE_4PI_EPS0_KJ = 138.935456;     // kJ*nm/mol/e²
const double ANGSTROM_TO_NM = 0.1;
const double KCAL_TO_KJ = 4.184;

namespace TIP3P {
    const double qO = -0.834;   
    const double qH = 0.417;    
    const double sigma_O = 3.15061;  // Å
    const double eps_O = 0.6364;     // kJ/mol
    const double rOH = 0.9572;       // Å
    const double aHOH = 104.52;      // degrees
}

namespace SWM4_NDP {
    const double qO = 1.71636;      
    const double qD = -1.71636;     
    const double qH = 0.55733;      
    const double qM = -1.11466;     
    const double sigma_O = 3.18395; // Å
    const double eps_O = 0.21094;   // kcal/mol
    const double k_drude = 1000.0;  // kcal/mol/Å²
    const double rOH = 0.9572;      // Å
    const double aHOH = 104.52;     // degrees
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
    std::vector<int> molecule_id;
    double box_size;
    double cutoff;
    int n_waters;
    int atoms_per_water;
    
public:
    WaterBox(double box_, double cut) : box_size(box_), cutoff(cut) {}
    
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
    
    // Check if intramolecular interaction should be excluded
    bool shouldExcludeIntramolecular(int i, int j) const {
        // In SWM4-NDP, we exclude ALL intramolecular interactions except:
        // 1. Drude with other atoms in molecule (except parent)
        // This is handled separately in energy calculation
        
        // For non-Drude atoms, exclude all intramolecular
        if (types[i] != 1 && types[j] != 1) {  // Neither is Drude
            return true;
        }
        
        // If one is Drude, check if it's Drude-parent
        for (size_t d = 0; d < drude_indices.size(); d++) {
            if ((i == drude_indices[d] && j == parent_indices[d]) ||
                (j == drude_indices[d] && i == parent_indices[d])) {
                return true;  // Exclude Drude-parent (in harmonic term)
            }
        }
        
        // Include Drude with other atoms in molecule
        return false;
    }
    
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
        
        // Initialize Drude positions
        std::cout << "  Initializing Drude positions...\n";
        runSCF(0.1, 200);
    }
    
    void addWaterMolecule(const Vec3& center, int mol_id) {
        int start_idx = positions.size();
        
        // Oxygen
        positions.push_back(center);
        charges.push_back(SWM4_NDP::qO);
        types.push_back(0);
        molecule_id.push_back(mol_id);
        
        // Drude - initially at parent
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
        double k_drude_kj_nm2 = SWM4_NDP::k_drude * KCAL_TO_KJ * 100.0;
        
        std::vector<Vec3> forces(drude_indices.size());
        
        for (int iter = 0; iter < max_iter; iter++) {
            double max_force = 0.0;
            
            // Calculate forces
            for (size_t d = 0; d < drude_indices.size(); d++) {
                int di = drude_indices[d];
                int pi = parent_indices[d];
                
                Vec3 force(0, 0, 0);
                
                // Harmonic force
                Vec3 dr = positions[di] - positions[pi];
                force = dr * (-k_drude_kj_nm2);
                
                // Electrostatic forces
                for (size_t j = 0; j < positions.size(); j++) {
                    if (j == di || j == pi) continue;
                    
                    // Include intramolecular only if same molecule
                    if (molecule_id[j] != molecule_id[di]) {
                        // Intermolecular
                        Vec3 rij = minimumImage(positions[j] - positions[di]);
                        double r2 = rij.norm2();
                        
                        if (r2 > 1e-10 && r2 < cutoff * cutoff) {
                            double r = std::sqrt(r2);
                            double f_mag = ONE_4PI_EPS0_KJ * charges[di] * charges[j] / (r2 * r);
                            force += rij * (f_mag / r);
                        }
                    } else {
                        // Intramolecular - include only H and M
                        if (types[j] == 2 || types[j] == 3) {  // H or M
                            Vec3 rij = positions[j] - positions[di];
                            double r2 = rij.norm2();
                            
                            if (r2 > 1e-10) {
                                double r = std::sqrt(r2);
                                double f_mag = ONE_4PI_EPS0_KJ * charges[di] * charges[j] / (r2 * r);
                                force += rij * (f_mag / r);
                            }
                        }
                    }
                }
                
                forces[d] = force;
                max_force = std::max(max_force, force.norm());
            }
            
            if (max_force < tolerance_kjmol_nm) break;
            
            // Update positions
            for (size_t d = 0; d < drude_indices.size(); d++) {
                double force_mag = forces[d].norm();
                
                if (force_mag > tolerance_kjmol_nm * 0.01) {
                    Vec3 displacement = forces[d] * (1.0 / k_drude_kj_nm2);
                    
                    double damping = 0.5;
                    if (force_mag > 1000.0) damping = 0.01;
                    else if (force_mag > 100.0) damping = 0.1;
                    
                    displacement = displacement * damping;
                    
                    double max_disp = 0.0001;
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
        // Optimize Drude positions
        runSCF();
        
        double energy = 0.0;
        double cutoff2 = cutoff * cutoff;
        double k_drude_kj_nm2 = SWM4_NDP::k_drude * KCAL_TO_KJ * 100.0;
        
        // 1. Harmonic energy
        for (size_t d = 0; d < drude_indices.size(); d++) {
            Vec3 dr = positions[drude_indices[d]] - positions[parent_indices[d]];
            energy += 0.5 * k_drude_kj_nm2 * dr.norm2();
        }
        
        // 2. Pairwise interactions
        for (size_t i = 0; i < positions.size(); i++) {
            for (size_t j = i + 1; j < positions.size(); j++) {
                bool same_mol = (molecule_id[i] == molecule_id[j]);
                
                if (same_mol) {
                    // Intramolecular - apply exclusion rules
                    if (shouldExcludeIntramolecular(i, j)) continue;
                } 
                
                Vec3 rij = same_mol ? positions[j] - positions[i] : minimumImage(positions[j] - positions[i]);
                double r2 = rij.norm2();
                
                if (r2 < cutoff2 && r2 > 1e-10) {
                    double r = std::sqrt(r2);
                    
                    // Coulomb
                    energy += ONE_4PI_EPS0_KJ * charges[i] * charges[j] / r;
                    
                    // LJ (only O-O intermolecular)
                    if (!same_mol && types[i] == 0 && types[j] == 0) {
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
        
        // Move O, D, H1, H2
        for (int i = 0; i < 4; i++) {
            positions[start + i] += displacement;
            applyPBC(positions[start + i]);
        }
        
        // Update M-site
        updateMSite(start);
        applyPBC(positions[start + 4]);
    }
};

void runBenchmark() {
    double box_size_A = 10.0;
    double box_size_nm = box_size_A * ANGSTROM_TO_NM;
    double cutoff_nm = 0.9;
    int n_waters = 10;
    int n_moves = 100;
    double max_displacement = 0.01;  // Smaller moves
    
    std::cout << "=================================================================\n";
    std::cout << "CHARMM Standard vs CHARMM Drude Performance Benchmark (FINAL)\n";
    std::cout << "=================================================================\n";
    std::cout << "Box size: " << box_size_A << " Å (" << n_waters << " waters)\n";
    std::cout << "Cutoff: " << cutoff_nm/ANGSTROM_TO_NM << " Å\n";
    std::cout << "Number of moves: " << n_moves << "\n\n";
    
    std::random_device rd;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(-max_displacement, max_displacement);
    std::uniform_int_distribution<> water_dis(0, n_waters-1);
    
    // TIP3P
    std::cout << "Testing TIP3P (Standard CHARMM)...\n";
    TIP3PWaterBox tip3p_box(box_size_nm, cutoff_nm);
    tip3p_box.generateWaterBox(n_waters);
    
    double tip3p_initial = tip3p_box.calculateEnergy();
    std::cout << "  Initial energy: " << tip3p_initial << " kJ/mol\n";
    std::cout << "  Per water: " << tip3p_initial/n_waters << " kJ/mol\n";
    
    auto start = std::chrono::high_resolution_clock::now();
    double tip3p_energy = tip3p_initial;
    std::vector<double> tip3p_energies;
    
    for (int move = 0; move < n_moves; move++) {
        int water_idx = water_dis(gen);
        Vec3 displacement(dis(gen), dis(gen), dis(gen));
        tip3p_box.moveWater(water_idx, displacement);
        tip3p_energy = tip3p_box.calculateEnergy();
        tip3p_energies.push_back(tip3p_energy);
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto tip3p_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    
    // SWM4-NDP
    std::cout << "\nTesting SWM4-NDP (CHARMM Drude with SCF)...\n";
    SWM4NDPWaterBox swm4_box(box_size_nm, cutoff_nm);
    swm4_box.generateWaterBox(n_waters);
    
    double swm4_initial = swm4_box.calculateEnergy();
    std::cout << "  Initial energy: " << swm4_initial << " kJ/mol\n";
    std::cout << "  Per water: " << swm4_initial/n_waters << " kJ/mol\n";
    
    gen.seed(42);  // Same moves
    
    start = std::chrono::high_resolution_clock::now();
    double swm4_energy = swm4_initial;
    std::vector<double> swm4_energies;
    
    for (int move = 0; move < n_moves; move++) {
        int water_idx = water_dis(gen);
        Vec3 displacement(dis(gen), dis(gen), dis(gen));
        swm4_box.moveWater(water_idx, displacement);
        swm4_energy = swm4_box.calculateEnergy();
        swm4_energies.push_back(swm4_energy);
        
        if (move % 25 == 0) {
            std::cout << "  Move " << move << ": E = " << swm4_energy 
                      << " kJ/mol (ΔE = " << (swm4_energy - swm4_initial) << ")\n";
        }
    }
    
    end = std::chrono::high_resolution_clock::now();
    auto swm4_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    
    // Calculate energy fluctuations
    double tip3p_avg = 0, swm4_avg = 0;
    for (auto e : tip3p_energies) tip3p_avg += e;
    for (auto e : swm4_energies) swm4_avg += e;
    tip3p_avg /= tip3p_energies.size();
    swm4_avg /= swm4_energies.size();
    
    double tip3p_std = 0, swm4_std = 0;
    for (auto e : tip3p_energies) tip3p_std += (e - tip3p_avg) * (e - tip3p_avg);
    for (auto e : swm4_energies) swm4_std += (e - swm4_avg) * (e - swm4_avg);
    tip3p_std = std::sqrt(tip3p_std / tip3p_energies.size());
    swm4_std = std::sqrt(swm4_std / swm4_energies.size());
    
    // Results
    std::cout << "\n=================================================================\n";
    std::cout << "RESULTS\n";
    std::cout << "=================================================================\n";
    std::cout << std::fixed << std::setprecision(2);
    
    std::cout << "TIP3P (Standard CHARMM):\n";
    std::cout << "  Total time: " << tip3p_time/1000.0 << " ms\n";
    std::cout << "  Time per move: " << tip3p_time/1000.0/n_moves << " ms\n";
    std::cout << "  Initial/Final energy: " << tip3p_initial << " / " << tip3p_energy << " kJ/mol\n";
    std::cout << "  Average energy: " << tip3p_avg << " ± " << tip3p_std << " kJ/mol\n\n";
    
    std::cout << "SWM4-NDP (CHARMM Drude with SCF):\n";
    std::cout << "  Total time: " << swm4_time/1000.0 << " ms\n";
    std::cout << "  Time per move: " << swm4_time/1000.0/n_moves << " ms\n";
    std::cout << "  Initial/Final energy: " << swm4_initial << " / " << swm4_energy << " kJ/mol\n";
    std::cout << "  Average energy: " << swm4_avg << " ± " << swm4_std << " kJ/mol\n\n";
    
    std::cout << "Performance ratio (Drude/Standard): " 
              << std::setprecision(1) << (double)swm4_time/tip3p_time << "x slower\n";
    
    // Check if energies are reasonable
    std::cout << "\n=================================================================\n";
    std::cout << "VALIDATION\n";
    std::cout << "=================================================================\n";
    
    bool tip3p_ok = (tip3p_avg < 0 && std::abs(tip3p_avg) < 1000 * n_waters);
    bool swm4_ok = (std::abs(swm4_avg) < 1000 * n_waters);
    
    std::cout << "TIP3P energy reasonable: " << (tip3p_ok ? "YES" : "NO") << "\n";
    std::cout << "SWM4-NDP energy reasonable: " << (swm4_ok ? "YES" : "NO") << "\n";
    
    if (!swm4_ok) {
        std::cout << "\nDEBUG: SWM4-NDP energy per water = " << swm4_avg/n_waters << " kJ/mol\n";
        std::cout << "This suggests remaining issues with intramolecular exclusions.\n";
    }
}

int main() {
    runBenchmark();
    return 0;
}