// Fixed version of Drude benchmark with proper SCF implementation

#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <random>
#include <iomanip>

const double ONE_4PI_EPS0 = 138.935456;  // kJ/mol·nm·e^-2
const double ANGSTROM_TO_NM = 0.1;
const double KCAL_TO_KJ = 4.184;

namespace SWM4_NDP {
    const double qO = 1.71636;    
    const double qD = -1.71636;   
    const double qH = 0.55733;    
    const double qM = -1.11466;   
    const double sigma_O = 0.318395;  // nm (already converted)
    const double eps_O = 0.21094 * KCAL_TO_KJ;  // kJ/mol
    const double k_drude = 1000.0 * KCAL_TO_KJ * 100.0;  // kcal/mol/Å² to kJ/mol/nm²
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

class SWM4NDPWaterBox {
private:
    std::vector<Vec3> positions;
    std::vector<double> charges;
    std::vector<int> types;  // 0=O, 1=D, 2=H, 3=M
    std::vector<int> drude_indices;
    std::vector<int> parent_indices;
    Vec3 box;
    double cutoff;
    int n_waters;
    const int atoms_per_water = 5;
    
public:
    SWM4NDPWaterBox(double box_size, double cutoff_) : 
        box(box_size, box_size, box_size), cutoff(cutoff_) {}
    
    void generateWaterBox(int target_waters) {
        n_waters = target_waters;
        
        // Clear existing data
        positions.clear();
        charges.clear();
        types.clear();
        drude_indices.clear();
        parent_indices.clear();
        
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
        
        // Initialize Drude positions properly
        initializeDrudePositions();
    }
    
    void addWaterMolecule(const Vec3& center) {
        int start_idx = positions.size();
        
        // Oxygen
        positions.push_back(center);
        charges.push_back(SWM4_NDP::qO);
        types.push_back(0);
        
        // Drude - initially at parent position
        positions.push_back(center);  // Will be optimized later
        charges.push_back(SWM4_NDP::qD);
        types.push_back(1);
        drude_indices.push_back(start_idx + 1);
        parent_indices.push_back(start_idx);
        
        // Hydrogen 1
        Vec3 h1_pos = center;
        h1_pos.x += SWM4_NDP::rOH;
        positions.push_back(h1_pos);
        charges.push_back(SWM4_NDP::qH);
        types.push_back(2);
        
        // Hydrogen 2  
        Vec3 h2_pos = center;
        h2_pos.x += SWM4_NDP::rOH * std::cos(SWM4_NDP::aHOH);
        h2_pos.y += SWM4_NDP::rOH * std::sin(SWM4_NDP::aHOH);
        positions.push_back(h2_pos);
        charges.push_back(SWM4_NDP::qH);
        types.push_back(2);
        
        // Virtual site M
        double w_o = 0.786646558;
        double w_h = 0.106676721;
        Vec3 m_pos;
        m_pos.x = w_o * center.x + w_h * h1_pos.x + w_h * h2_pos.x;
        m_pos.y = w_o * center.y + w_h * h1_pos.y + w_h * h2_pos.y;
        m_pos.z = w_o * center.z + w_h * h1_pos.z + w_h * h2_pos.z;
        positions.push_back(m_pos);
        charges.push_back(SWM4_NDP::qM);
        types.push_back(3);
    }
    
    void initializeDrudePositions() {
        // Two-stage initialization for better stability
        
        // Stage 1: Each Drude in isolation (no intermolecular forces)
        for (size_t d = 0; d < drude_indices.size(); d++) {
            int di = drude_indices[d];
            int pi = parent_indices[d];
            
            // For isolated Drude, equilibrium is at parent
            positions[di] = positions[pi];
        }
        
        // Stage 2: Include intermolecular interactions
        runSCF(1.0, 200);  // More iterations for initial setup
    }
    
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
    
    void runSCF(double tolerance = 1.0, int max_iter = 50) {
        std::vector<Vec3> forces(drude_indices.size());
        
        for (int iter = 0; iter < max_iter; iter++) {
            double max_force = 0.0;
            
            // Calculate forces on all Drudes
            for (size_t d = 0; d < drude_indices.size(); d++) {
                int di = drude_indices[d];
                int pi = parent_indices[d];
                
                Vec3 force(0, 0, 0);
                
                // Harmonic restoring force
                Vec3 dr = positions[di] - positions[pi];
                force = force - (dr * SWM4_NDP::k_drude);
                
                // Coulomb forces from all atoms
                for (size_t j = 0; j < positions.size(); j++) {
                    if (j == di) continue;  // Skip self
                    
                    // Skip parent (force already included in harmonic term)
                    if (j == pi) continue;
                    
                    Vec3 rij = minimumImage(positions[j] - positions[di]);
                    double r2 = rij.norm2();
                    
                    if (r2 < cutoff * cutoff && r2 > 1e-10) {
                        double r = std::sqrt(r2);
                        double f_mag = ONE_4PI_EPS0 * charges[di] * charges[j] / (r2 * r);
                        force = force - (rij * (f_mag / r));
                    }
                }
                
                forces[d] = force;
                max_force = std::max(max_force, force.norm());
            }
            
            if (max_force < tolerance) {
                break;
            }
            
            // Update positions with strict limits
            for (size_t d = 0; d < drude_indices.size(); d++) {
                double force_mag = forces[d].norm();
                
                if (force_mag > tolerance * 0.01) {
                    // Adaptive damping
                    double damping;
                    if (force_mag > 1000.0) {
                        damping = 0.001;
                    } else if (force_mag > 100.0) {
                        damping = 0.01;
                    } else if (force_mag > 10.0) {
                        damping = 0.1;
                    } else {
                        damping = 0.3;
                    }
                    
                    Vec3 delta = forces[d] * (damping / SWM4_NDP::k_drude);
                    
                    // Hard limit on displacement
                    const double max_disp = 0.0001;  // 0.001 Å
                    double delta_mag = delta.norm();
                    if (delta_mag > max_disp) {
                        delta = delta * (max_disp / delta_mag);
                    }
                    
                    positions[drude_indices[d]] = positions[drude_indices[d]] + delta;
                }
            }
        }
    }
    
    double calculateEnergy() {
        // First optimize Drude positions
        runSCF();
        
        double energy = 0.0;
        double cutoff2 = cutoff * cutoff;
        
        // Harmonic energy for Drude oscillators
        for (size_t d = 0; d < drude_indices.size(); d++) {
            Vec3 dr = positions[drude_indices[d]] - positions[parent_indices[d]];
            energy += 0.5 * SWM4_NDP::k_drude * dr.norm2();
        }
        
        // Pairwise interactions
        for (size_t i = 0; i < positions.size(); i++) {
            for (size_t j = i + 1; j < positions.size(); j++) {
                // Skip intramolecular interactions
                if (i/atoms_per_water == j/atoms_per_water) continue;
                
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
                
                if (r2 < cutoff2 && r2 > 1e-10) {
                    double r = std::sqrt(r2);
                    
                    // Coulomb
                    energy += ONE_4PI_EPS0 * charges[i] * charges[j] / r;
                    
                    // LJ (only O-O)
                    if (types[i] == 0 && types[j] == 0) {
                        double r6 = r2 * r2 * r2;
                        double sigma6 = SWM4_NDP::sigma_O * SWM4_NDP::sigma_O * SWM4_NDP::sigma_O;
                        sigma6 = sigma6 * sigma6;
                        double r12 = r6 * r6;
                        energy += 4.0 * SWM4_NDP::eps_O * (sigma6*sigma6/r12 - sigma6/r6);
                    }
                }
            }
        }
        
        return energy;
    }
    
    void moveWater(int water_idx, const Vec3& displacement) {
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
    }
};

int main() {
    // Test with a small system first
    double box_size = 10.0 * ANGSTROM_TO_NM;  // 10 Å
    double cutoff = 0.9;  // 9 Å
    int n_waters = 3;  // Just 3 waters for testing
    
    std::cout << "=== Testing Fixed SWM4-NDP Implementation ===\n";
    std::cout << "Box size: " << box_size/ANGSTROM_TO_NM << " Å\n";
    std::cout << "Number of waters: " << n_waters << "\n";
    std::cout << "k_drude = " << SWM4_NDP::k_drude << " kJ/mol/nm²\n\n";
    
    SWM4NDPWaterBox water_box(box_size, cutoff);
    water_box.generateWaterBox(n_waters);
    
    std::cout << "Calculating initial energy...\n";
    double initial_energy = water_box.calculateEnergy();
    std::cout << "Initial energy: " << initial_energy << " kJ/mol\n";
    std::cout << "Energy per water: " << initial_energy/n_waters << " kJ/mol\n\n";
    
    // Test a few moves
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.01, 0.01);  // Small moves
    
    std::cout << "Testing 10 MC moves...\n";
    for (int move = 0; move < 10; move++) {
        int water_idx = move % n_waters;
        Vec3 displacement(dis(gen), dis(gen), dis(gen));
        water_box.moveWater(water_idx, displacement);
        
        double energy = water_box.calculateEnergy();
        std::cout << "Move " << move+1 << ": E = " << energy << " kJ/mol\n";
    }
    
    return 0;
}