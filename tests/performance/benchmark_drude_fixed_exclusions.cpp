// Fixed SWM4-NDP implementation with correct intramolecular exclusions

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
    
    // Get atom type within molecule (0=O, 1=D, 2=H, 3=M)
    int getAtomTypeInMolecule(int global_idx) const {
        int mol_id = molecule_id[global_idx];
        int mol_start = mol_id * atoms_per_water;
        return global_idx - mol_start;
    }
    
    // Check if intramolecular interaction should be included
    bool shouldIncludeIntramolecular(int i, int j) const {
        if (molecule_id[i] != molecule_id[j]) return true;  // Different molecules
        
        int type_i = getAtomTypeInMolecule(i);
        int type_j = getAtomTypeInMolecule(j);
        
        // Sort to make checking easier
        if (type_i > type_j) std::swap(type_i, type_j);
        
        // EXCLUSIONS (return false for these):
        // O-D (0-1): Excluded (in harmonic term)
        // O-H (0-2): Excluded (bonded)
        // O-M (0-3): Excluded (M is virtual site of O,H,H)
        // D-D: Not possible (only one D per molecule)
        // H-H (2-2): Excluded (1-3 connected)
        // H-M (2-3): Excluded (M is virtual site)
        
        // INCLUSIONS (return true for these):
        // D-H (1-2): Include
        // D-M (1-3): Include
        
        // Check specific pairs
        if (type_i == 0 && type_j == 1) return false;  // O-D
        if (type_i == 0 && type_j == 2) return false;  // O-H
        if (type_i == 0 && type_j == 3) return false;  // O-M
        if (type_i == 2 && type_j == 2) return false;  // H-H
        if (type_i == 2 && type_j == 3) return false;  // H-M
        
        if (type_i == 1 && type_j == 2) return true;   // D-H
        if (type_i == 1 && type_j == 3) return true;   // D-M
        
        return false;  // Default: exclude
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
        
        // Debug: Check one water molecule
        std::cout << "  First water molecule charges and positions:\n";
        for (int i = 0; i < 5; i++) {
            std::cout << "    Atom " << i << " (type " << types[i] << "): "
                      << "q=" << charges[i] << ", pos=(" 
                      << positions[i].x << ", " << positions[i].y << ", " << positions[i].z << ")\n";
        }
    }
    
    void addWaterMolecule(const Vec3& center, int mol_id) {
        int start_idx = positions.size();
        
        // Oxygen (type 0)
        positions.push_back(center);
        charges.push_back(SWM4_NDP::qO);
        types.push_back(0);
        molecule_id.push_back(mol_id);
        
        // Drude (type 1) - initially at parent
        positions.push_back(center);
        charges.push_back(SWM4_NDP::qD);
        types.push_back(1);
        molecule_id.push_back(mol_id);
        drude_indices.push_back(start_idx + 1);
        parent_indices.push_back(start_idx);
        
        // Hydrogen 1 (type 2)
        double angle_rad = SWM4_NDP::aHOH * M_PI / 180.0;
        Vec3 h1 = center;
        h1.x += SWM4_NDP::rOH * ANGSTROM_TO_NM;
        positions.push_back(h1);
        charges.push_back(SWM4_NDP::qH);
        types.push_back(2);
        molecule_id.push_back(mol_id);
        
        // Hydrogen 2 (type 2)
        Vec3 h2 = center;
        h2.x += SWM4_NDP::rOH * ANGSTROM_TO_NM * std::cos(angle_rad);
        h2.y += SWM4_NDP::rOH * ANGSTROM_TO_NM * std::sin(angle_rad);
        positions.push_back(h2);
        charges.push_back(SWM4_NDP::qH);
        types.push_back(2);
        molecule_id.push_back(mol_id);
        
        // Virtual M-site (type 3)
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
            
            // Calculate forces on all Drudes
            for (size_t d = 0; d < drude_indices.size(); d++) {
                int di = drude_indices[d];
                int pi = parent_indices[d];
                
                Vec3 force(0, 0, 0);
                
                // Harmonic force
                Vec3 dr = positions[di] - positions[pi];
                force = dr * (-k_drude_kj_nm2);
                
                // Electrostatic forces
                for (size_t j = 0; j < positions.size(); j++) {
                    if (j == di) continue;  // Skip self
                    
                    // Check if this interaction should be included
                    if (!shouldIncludeIntramolecular(di, j)) continue;
                    
                    Vec3 rij;
                    if (molecule_id[di] == molecule_id[j]) {
                        rij = positions[j] - positions[di];  // No PBC for intramolecular
                    } else {
                        rij = minimumImage(positions[j] - positions[di]);  // PBC for intermolecular
                    }
                    
                    double r2 = rij.norm2();
                    
                    if (r2 > 1e-10 && (molecule_id[di] == molecule_id[j] || r2 < cutoff * cutoff)) {
                        double r = std::sqrt(r2);
                        double f_mag = ONE_4PI_EPS0_KJ * charges[di] * charges[j] / (r2 * r);
                        force += rij * (f_mag / r);
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
                // Check if this interaction should be included
                if (!shouldIncludeIntramolecular(i, j)) continue;
                
                Vec3 rij;
                if (molecule_id[i] == molecule_id[j]) {
                    rij = positions[j] - positions[i];  // No PBC for intramolecular
                } else {
                    rij = minimumImage(positions[j] - positions[i]);  // PBC for intermolecular
                }
                
                double r2 = rij.norm2();
                
                // For intramolecular, no cutoff; for intermolecular, apply cutoff
                bool should_calculate = (molecule_id[i] == molecule_id[j]) || (r2 < cutoff2);
                
                if (should_calculate && r2 > 1e-10) {
                    double r = std::sqrt(r2);
                    
                    // Coulomb
                    energy += ONE_4PI_EPS0_KJ * charges[i] * charges[j] / r;
                    
                    // LJ (only O-O intermolecular)
                    if (molecule_id[i] != molecule_id[j] && types[i] == 0 && types[j] == 0) {
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

int main() {
    double box_size_A = 10.0;
    double box_size_nm = box_size_A * ANGSTROM_TO_NM;
    double cutoff_nm = 0.9;
    int n_waters = 3;  // Start with just 3 waters
    
    std::cout << "=================================================================\n";
    std::cout << "Testing SWM4-NDP with Correct Exclusions\n";
    std::cout << "=================================================================\n";
    std::cout << "Box size: " << box_size_A << " Å (" << n_waters << " waters)\n";
    std::cout << "Cutoff: " << cutoff_nm/ANGSTROM_TO_NM << " Å\n\n";
    
    // First test single water intramolecular energy
    std::cout << "Single water molecule test:\n";
    {
        SWM4NDPWaterBox single_water(100.0, 50.0);  // Large box
        single_water.generateWaterBox(1);
        double single_energy = single_water.calculateEnergy();
        std::cout << "  Single water energy: " << single_energy << " kJ/mol\n";
        std::cout << "  (Should be small, mostly from D-H and D-M interactions)\n\n";
    }
    
    // Test with 3 waters
    std::cout << "Three water molecules test:\n";
    SWM4NDPWaterBox water_box(box_size_nm, cutoff_nm);
    water_box.generateWaterBox(n_waters);
    
    double initial_energy = water_box.calculateEnergy();
    std::cout << "  Initial energy: " << initial_energy << " kJ/mol\n";
    std::cout << "  Per water: " << initial_energy/n_waters << " kJ/mol\n";
    
    // Test a few moves
    std::random_device rd;
    std::mt19937 gen(42);
    std::uniform_real_distribution<> dis(-0.01, 0.01);
    
    std::cout << "\nTesting 10 MC moves:\n";
    for (int move = 0; move < 10; move++) {
        int water_idx = move % n_waters;
        Vec3 displacement(dis(gen), dis(gen), dis(gen));
        water_box.moveWater(water_idx, displacement);
        
        double energy = water_box.calculateEnergy();
        std::cout << "  Move " << move+1 << ": E = " << energy 
                  << " kJ/mol (ΔE = " << (energy - initial_energy) << ")\n";
    }
    
    return 0;
}