// examples/basic/cpp_example.cpp

#include "pygcmc/core/system.hpp"
#include <iostream>

int main(int argc, char** argv) {
    if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <pdb_file> <psf_file> <top_file> <forcefield_file>" << std::endl;
        return 1;
    }

    std::string pdb_file = argv[1];
    std::string psf_file = argv[2];
    std::string top_file = argv[3];
    std::string ff_file = argv[4];

    pygcmc::core::System system;

    // 加载文件
    system.load_pdb(pdb_file);
    system.load_psf(psf_file);
    system.load_top(top_file);
    system.load_forcefield(ff_file);

    // 计算总能量
    double energy = system.compute_total_energy();
    std::cout << "Total Energy: " << energy << std::endl;

    return 0;
}
