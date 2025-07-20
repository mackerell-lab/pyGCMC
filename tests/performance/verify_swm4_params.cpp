// Verify SWM4-NDP parameters from literature

#include <iostream>
#include <cmath>
#include <iomanip>

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "=== SWM4-NDP Parameter Verification ===\n\n";
    
    // From Lamoureux et al. 2003 and other sources
    std::cout << "Literature values for SWM4-NDP water:\n";
    std::cout << "  Oxygen charge: +1.71636 e\n";
    std::cout << "  Drude charge: -1.71636 e\n";
    std::cout << "  Hydrogen charge: +0.55733 e\n";
    std::cout << "  M-site charge: -1.11466 e\n";
    std::cout << "  Total charge: 0.0 e\n\n";
    
    std::cout << "Force field parameters:\n";
    std::cout << "  k_D = 1000 kcal/mol/Å² (Drude force constant)\n";
    std::cout << "  α_O = 0.978 Å³ (oxygen polarizability)\n";
    std::cout << "  σ_O = 3.18395 Å (LJ sigma for oxygen)\n";
    std::cout << "  ε_O = 0.21094 kcal/mol (LJ epsilon for oxygen)\n\n";
    
    // Calculate what the force constant should be
    const double alpha_A3 = 0.978;  // Å³
    const double qD = -1.71636;
    const double e2_per_A = 332.0637;  // e²/Å in kcal/mol
    
    std::cout << "Force constant from polarizability:\n";
    std::cout << "  k = q²/(4πε₀α)\n";
    std::cout << "  Using e²/(4πε₀) = " << e2_per_A << " kcal·Å/mol\n";
    std::cout << "  k = " << qD*qD << " × " << e2_per_A << " / " << alpha_A3 << "\n";
    std::cout << "  k = " << (qD*qD * e2_per_A / alpha_A3) << " kcal/mol/Å²\n\n";
    
    // The issue: we should NOT use the force constant directly!
    // Instead, we should use the polarizability
    
    std::cout << "IMPORTANT REALIZATION:\n";
    std::cout << "In the Drude model, we specify EITHER:\n";
    std::cout << "1. The force constant k, OR\n";
    std::cout << "2. The polarizability α\n";
    std::cout << "They are related by: k = q²/(4πε₀α)\n\n";
    
    std::cout << "For SWM4-NDP, the standard approach is:\n";
    std::cout << "- Specify α = 0.978 Å³\n";
    std::cout << "- Calculate k from α and q\n";
    std::cout << "- This gives k ≈ 1000 kcal/mol/Å²\n\n";
    
    // Check energy scale
    double r_test = 0.01;  // Å displacement
    double k_kcal = 1000.0;
    double E_kcal = 0.5 * k_kcal * r_test * r_test;
    std::cout << "Energy for " << r_test << " Å displacement:\n";
    std::cout << "  E = 0.5 × " << k_kcal << " × " << r_test << "² = " << E_kcal << " kcal/mol\n";
    std::cout << "  E = " << E_kcal * 4.184 << " kJ/mol\n";
    
    return 0;
}