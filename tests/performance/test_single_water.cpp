// Test single SWM4-NDP water to verify basic physics

#include <iostream>
#include <cmath>
#include <iomanip>

const double ONE_4PI_EPS0 = 138.935456;
const double KCAL_TO_KJ = 4.184;

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    // SWM4-NDP parameters
    const double qD = -1.71636;
    const double k_drude_kcal = 1000.0;  // kcal/mol/Å²
    
    std::cout << "=== SWM4-NDP Single Drude Analysis ===\n\n";
    
    // Convert force constant
    std::cout << "Force constant conversions:\n";
    std::cout << "  Original: " << k_drude_kcal << " kcal/mol/Å²\n";
    
    // Method 1: Direct conversion
    double k_drude_kj_nm2 = k_drude_kcal * KCAL_TO_KJ * 100.0;  // 1 Å = 0.1 nm, so 1/Å² = 100/nm²
    std::cout << "  In kJ/mol/nm²: " << k_drude_kj_nm2 << "\n";
    
    // Calculate expected equilibrium distance for isolated Drude
    // At equilibrium: F_harmonic = 0 (no external field)
    // So Drude should be at parent position
    
    std::cout << "\nFor isolated Drude (no external field):\n";
    std::cout << "  Expected position: at parent (r = 0)\n";
    
    // Calculate polarizability
    double alpha = qD * qD / (ONE_4PI_EPS0 * k_drude_kj_nm2);
    std::cout << "\nPolarizability:\n";
    std::cout << "  α = q²/(4πε₀k) = " << alpha << " nm³\n";
    std::cout << "  α = " << alpha * 1000.0 << " Å³\n";
    
    // Test energy for various displacements
    std::cout << "\nHarmonic energy for various displacements:\n";
    std::cout << "  r (Å)    r (nm)     E (kJ/mol)\n";
    std::cout << "  ------   --------   -----------\n";
    
    double r_values[] = {0.0001, 0.001, 0.01, 0.1, 1.0};  // Å
    for (double r_ang : r_values) {
        double r_nm = r_ang * 0.1;
        double E = 0.5 * k_drude_kj_nm2 * r_nm * r_nm;
        std::cout << "  " << std::setw(6) << r_ang 
                  << "   " << std::setw(8) << r_nm 
                  << "   " << std::setw(11) << E << "\n";
    }
    
    // Expected behavior with external field
    std::cout << "\nWith external electric field:\n";
    std::cout << "  Induced dipole μ = αE\n";
    std::cout << "  For E = 1 V/nm, μ = " << alpha << " e·nm\n";
    std::cout << "  Drude displacement = μ/|q| = " << alpha/std::abs(qD) << " nm\n";
    
    return 0;
}