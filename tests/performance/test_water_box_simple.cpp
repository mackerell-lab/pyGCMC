#include <iostream>
#include <chrono>
#include <random>
#include <vector>
#include <cmath>
#include <fstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Include the actual headers used in bindings
#include "model/montecarlo/MCStructures.hpp"
#include "simulation/simulation.hpp"

namespace py = pybind11;
using namespace pygcmc;
using namespace std;
using namespace std::chrono;

int main() {
    cout << "============================================================\n";
    cout << "Water Box Performance Test - C++ Version (Non-Drude)\n";
    cout << "============================================================\n";
    
    // Initialize Python interpreter for pybind11
    py::scoped_interpreter guard{};
    
    try {
        // Parameters
        const double boxSize = 3.0;  // nm
        const double waterSpacing = 0.31;  // nm
        const int nSteps = 20000;
        const double maxDisplacement = 0.001;  // nm per step
        const double maxRotation = 0.01;  // radians per step
        
        // Create state
        model::MCState state;
        state.info.box[0] = boxSize;
        state.info.box[1] = boxSize;
        state.info.box[2] = boxSize;
        state.info.cutoff = min(boxSize/2 - 0.1, 1.2);
        state.info.temperature = 300.0;
        
        // Calculate number of waters
        int nPerDim = int(boxSize / waterSpacing);
        int totalWaters = nPerDim * nPerDim * nPerDim;
        
        cout << "Creating water box:\n";
        cout << "  Box size: " << boxSize << " nm\n";
        cout << "  Waters per dimension: " << nPerDim << "\n";
        cout << "  Total water molecules: " << totalWaters << "\n";
        cout << "  Total atoms: " << totalWaters * 3 << "\n";
        cout << "  Cutoff: " << state.info.cutoff << " nm\n";
        
        // Random number generators
        random_device rd;
        mt19937 gen(42);  // Fixed seed
        uniform_real_distribution<> smallDisp(-0.01, 0.01);
        uniform_real_distribution<> dispDist(-maxDisplacement, maxDisplacement);
        uniform_real_distribution<> rotDist(-maxRotation, maxRotation);
        
        // Create water molecules
        int molId = 0;
        for (int ix = 0; ix < nPerDim; ix++) {
            for (int iy = 0; iy < nPerDim; iy++) {
                for (int iz = 0; iz < nPerDim; iz++) {
                    double x = (ix + 0.5) * waterSpacing + smallDisp(gen);
                    double y = (iy + 0.5) * waterSpacing + smallDisp(gen);
                    double z = (iz + 0.5) * waterSpacing + smallDisp(gen);
                    
                    // Create water - oxygen
                    model::MCAtom oxygen;
                    oxygen.x = x;
                    oxygen.y = y;
                    oxygen.z = z;
                    oxygen.charge = -0.834;  // TIP3P
                    oxygen.type = 0;
                    state.atoms.push_back(oxygen);
                    
                    // Hydrogen 1
                    model::MCAtom h1;
                    h1.x = x + 0.09572;
                    h1.y = y;
                    h1.z = z;
                    h1.charge = 0.417;
                    h1.type = 1;
                    state.atoms.push_back(h1);
                    
                    // Hydrogen 2
                    model::MCAtom h2;
                    double angle = 104.52 * M_PI / 180.0;
                    h2.x = x + 0.09572 * cos(angle);
                    h2.y = y + 0.09572 * sin(angle);
                    h2.z = z;
                    h2.charge = 0.417;
                    h2.type = 1;
                    state.atoms.push_back(h2);
                    
                    // Create residue
                    model::MCResidue res;
                    res.atomStart = molId * 3;
                    res.atomCount = 3;
                    res.active = true;
                    res.fixed = false;
                    res.type = 0;
                    state.residues.push_back(res);
                    
                    molId++;
                }
            }
        }
        
        state.activeAtomCount = state.atoms.size();
        state.activeResidueCount = state.residues.size();
        
        // Set force field - TIP3P
        state.forcefield.numTotalTypes = 2;
        state.forcefield.numMovementTypes = 2;
        
        // LJ parameters
        state.forcefield.ljSigma = {0.315, 0.0, 0.0, 0.0};
        state.forcefield.ljEps = {0.636, 0.0, 0.0, 0.0};
        
        // Select molecule to move (center one)
        int targetMol = totalWaters / 2;
        int startAtom = state.residues[targetMol].atomStart;
        
        cout << "\nTesting molecule " << targetMol << " movement\n";
        cout << "Initial O position: (" 
             << state.atoms[startAtom].x << ", " 
             << state.atoms[startAtom].y << ", " 
             << state.atoms[startAtom].z << ")\n";
        
        // Calculate initial energy
        auto start = high_resolution_clock::now();
        simulation::computeSystemEnergyCutoff(state);
        auto end = high_resolution_clock::now();
        
        double initialEnergy = 0.0;
        for (const auto& res : state.residues) {
            initialEnergy += res.energy_elec + res.energy_vdw;
        }
        
        auto initialTime = duration_cast<microseconds>(end - start).count() / 1000.0;
        cout << "Initial energy: " << initialEnergy << " kJ/mol\n";
        cout << "Initial calculation time: " << initialTime << " ms\n";
        
        // Arrays for statistics
        vector<double> energies;
        vector<double> times;
        energies.reserve(nSteps);
        times.reserve(nSteps);
        
        // Progress markers
        vector<int> progressSteps = {2000, 5000, 10000, 15000, 18000, 20000};
        
        cout << "\nRunning " << nSteps << " steps...\n";
        auto totalStart = high_resolution_clock::now();
        
        // Main simulation loop
        for (int step = 0; step < nSteps; step++) {
            // Random translation
            double dx = dispDist(gen);
            double dy = dispDist(gen);
            double dz = dispDist(gen);
            
            // Random rotation angle
            double angle = rotDist(gen);
            double cosA = cos(angle);
            double sinA = sin(angle);
            
            // Calculate center of mass
            double comX = 0, comY = 0, comZ = 0;
            for (int i = 0; i < 3; i++) {
                comX += state.atoms[startAtom + i].x;
                comY += state.atoms[startAtom + i].y;
                comZ += state.atoms[startAtom + i].z;
            }
            comX /= 3.0;
            comY /= 3.0;
            comZ /= 3.0;
            
            // Move and rotate water
            for (int i = 0; i < 3; i++) {
                auto& atom = state.atoms[startAtom + i];
                
                // Translate to origin
                double relX = atom.x - comX;
                double relY = atom.y - comY;
                
                // Rotate around z-axis
                double newRelX = relX * cosA - relY * sinA;
                double newRelY = relX * sinA + relY * cosA;
                
                // Translate back and add displacement
                atom.x = newRelX + comX + dx;
                atom.y = newRelY + comY + dy;
                atom.z += dz;
                
                // Apply PBC
                atom.x = atom.x - floor(atom.x / boxSize) * boxSize;
                atom.y = atom.y - floor(atom.y / boxSize) * boxSize;
                atom.z = atom.z - floor(atom.z / boxSize) * boxSize;
            }
            
            // Calculate energy
            auto stepStart = high_resolution_clock::now();
            simulation::computeSystemEnergyCutoff(state);
            auto stepEnd = high_resolution_clock::now();
            
            double energy = 0.0;
            for (const auto& res : state.residues) {
                energy += res.energy_elec + res.energy_vdw;
            }
            
            double stepTime = duration_cast<microseconds>(stepEnd - stepStart).count() / 1000.0;
            energies.push_back(energy);
            times.push_back(stepTime);
            
            // Progress update
            if (find(progressSteps.begin(), progressSteps.end(), step + 1) != progressSteps.end()) {
                double avgTime = 0;
                for (double t : times) avgTime += t;
                avgTime /= times.size();
                
                cout << "  Progress: " << (step + 1) * 100 / nSteps << "% (" 
                     << (step + 1) << "/" << nSteps << " steps)\n";
                cout << "    Current energy: " << energy << " kJ/mol\n";
                cout << "    Avg time per step: " << avgTime << " ms\n";
            }
        }
        
        auto totalEnd = high_resolution_clock::now();
        double totalTime = duration_cast<milliseconds>(totalEnd - totalStart).count() / 1000.0;
        
        // Calculate statistics
        double avgEnergy = 0, avgTime = 0;
        for (size_t i = 0; i < energies.size(); i++) {
            avgEnergy += energies[i];
            avgTime += times[i];
        }
        avgEnergy /= energies.size();
        avgTime /= times.size();
        
        cout << "\n=== Performance Summary ===\n";
        cout << "Total simulation time: " << totalTime << " seconds\n";
        cout << "Average time per step: " << avgTime << " ms\n";
        cout << "Steps per second: " << nSteps / totalTime << "\n";
        
        cout << "\n=== Energy Statistics ===\n";
        cout << "Average energy: " << avgEnergy << " kJ/mol\n";
        cout << "Energy drift: " << energies.back() - energies.front() << " kJ/mol\n";
        
        // Save results
        ofstream outFile("water_box_performance_cpp.txt");
        outFile << "Box size: " << boxSize << " nm\n";
        outFile << "Total waters: " << totalWaters << "\n";
        outFile << "Total atoms: " << totalWaters * 3 << "\n";
        outFile << "Steps: " << nSteps << "\n";
        outFile << "Avg time per step: " << avgTime << " ms\n";
        outFile << "Total time: " << totalTime << " s\n";
        outFile << "Avg energy: " << avgEnergy << " kJ/mol\n";
        outFile.close();
        
        cout << "\nResults saved to water_box_performance_cpp.txt\n";
        
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}