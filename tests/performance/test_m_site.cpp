// Test M-site calculation

#include <iostream>
#include <cmath>
#include <iomanip>

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    // SWM4-NDP geometry (in Å)
    double rOH = 0.9572;
    double aHOH = 104.52 * M_PI / 180.0;  // radians
    
    // Positions in Å
    double O_x = 0.0, O_y = 0.0, O_z = 0.0;
    double H1_x = rOH, H1_y = 0.0, H1_z = 0.0;
    double H2_x = rOH * std::cos(aHOH), H2_y = rOH * std::sin(aHOH), H2_z = 0.0;
    
    std::cout << "=== SWM4-NDP M-site Calculation ===\n\n";
    std::cout << "Atomic positions (Å):\n";
    std::cout << "  O:  (" << O_x << ", " << O_y << ", " << O_z << ")\n";
    std::cout << "  H1: (" << H1_x << ", " << H1_y << ", " << H1_z << ")\n";
    std::cout << "  H2: (" << H2_x << ", " << H2_y << ", " << H2_z << ")\n\n";
    
    // M-site weights
    double w_O = 0.786646558;
    double w_H = 0.106676721;
    
    std::cout << "M-site weights:\n";
    std::cout << "  w_O = " << w_O << "\n";
    std::cout << "  w_H = " << w_H << " (each H)\n";
    std::cout << "  Sum = " << (w_O + 2*w_H) << " (should be 1.0)\n\n";
    
    // Calculate M position
    double M_x = w_O * O_x + w_H * H1_x + w_H * H2_x;
    double M_y = w_O * O_y + w_H * H1_y + w_H * H2_y;
    double M_z = w_O * O_z + w_H * H1_z + w_H * H2_z;
    
    std::cout << "M-site position:\n";
    std::cout << "  M: (" << M_x << ", " << M_y << ", " << M_z << ") Å\n\n";
    
    // Distances
    double r_OM = std::sqrt(M_x*M_x + M_y*M_y + M_z*M_z);
    double r_H1M = std::sqrt((H1_x-M_x)*(H1_x-M_x) + (H1_y-M_y)*(H1_y-M_y) + (H1_z-M_z)*(H1_z-M_z));
    double r_H2M = std::sqrt((H2_x-M_x)*(H2_x-M_x) + (H2_y-M_y)*(H2_y-M_y) + (H2_z-M_z)*(H2_z-M_z));
    
    std::cout << "Distances (Å):\n";
    std::cout << "  O-M:  " << r_OM << "\n";
    std::cout << "  H1-M: " << r_H1M << "\n";
    std::cout << "  H2-M: " << r_H2M << "\n\n";
    
    // Energy scale check
    double qO = 1.71636;
    double qM = -1.11466;
    double ONE_4PI_EPS0 = 332.0637;  // kcal*Å/mol
    
    double E_OM = ONE_4PI_EPS0 * qO * qM / r_OM;
    std::cout << "O-M Coulomb energy: " << E_OM << " kcal/mol\n";
    std::cout << "                   = " << E_OM * 4.184 << " kJ/mol\n\n";
    
    // Verify charge neutrality  
    double qH = 0.55733;
    double qD = -1.71636;
    double q_total = qO + qD + 2*qH + qM;
    std::cout << "Total charge: " << q_total << " (should be 0)\n";
    
    return 0;
}