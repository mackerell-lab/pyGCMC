#include "PGPContext.hpp"
#include "PGPCore.hpp"
#include "../pme/PMECore.hpp"
#include "../pme/PMEGridCharge.hpp"
#include "../pme/PMERecip.hpp"
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <complex>

namespace platform {
namespace cpu {
namespace energy {
namespace pgp {

void PGPContext::initialize(double cutoff, 
                          const std::array<double, 3>& box,
                          double alpha,
                          const std::array<int, 3>& meshSize,
                          double potential_cutoff,
                          const std::array<int, 3>& potentialGridSize,
                          int splineOrder,
                          double tolerance) {
    // Store basic parameters
    params_.cutoff = cutoff;
    params_.box = box;
    params_.alpha = alpha;
    params_.meshSize = meshSize;
    params_.potential_cutoff = potential_cutoff;
    params_.potentialGridSize = potentialGridSize;
    params_.splineOrder = splineOrder;
    params_.tolerance = tolerance;
    
    // Calculate derived parameters
    params_.alphaEwald = alpha;
    params_.ewaldCoeff = -params_.alphaEwald * params_.alphaEwald;
    params_.selfEnergyCoeff = -params_.alphaEwald / std::sqrt(M_PI) * 138.935456;
    
    // Setup grid dimensions
    params_.gridSizeX = meshSize[0];
    params_.gridSizeY = meshSize[1];
    params_.gridSizeZ = meshSize[2];
    params_.gridSizeYZ = params_.gridSizeY * params_.gridSizeZ;
    params_.gridTotal = params_.gridSizeX * params_.gridSizeYZ;
    
    // Initialize tables and grids
    setupTables();
    setupGrids();
}

void PGPContext::setupTables() {
    // Initialize ERFC table
    params_.erfcTable.resize(Parameters::ERFC_TABLE_SIZE);
    for (int i = 0; i < Parameters::ERFC_TABLE_SIZE; ++i) {
        double x = i / Parameters::ERFC_TABLE_SCALE;
        params_.erfcTable[i] = std::erfc(x);
    }
    
    // Initialize EXP table
    params_.expTable.resize(Parameters::EXP_TABLE_SIZE);
    params_.dExpTable.resize(Parameters::EXP_TABLE_SIZE);
    for (int i = 0; i < Parameters::EXP_TABLE_SIZE; ++i) {
        double x = -i * i / (Parameters::EXP_TABLE_SCALE * Parameters::EXP_TABLE_SCALE);
        params_.expTable[i] = std::exp(x);
        params_.dExpTable[i] = -2.0 * i / (Parameters::EXP_TABLE_SCALE * Parameters::EXP_TABLE_SCALE) * params_.expTable[i];
    }
}

void PGPContext::setupGrids() {
    // Allocate grids
    params_.pmeGrid.resize(params_.gridTotal);
    params_.pmeGridSaved.resize(params_.gridTotal);
    params_.potentialGrid.resize(params_.gridTotal * 10); // Assuming max 10 atom types
    
    // Initialize to zero
    std::fill(params_.pmeGrid.begin(), params_.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    std::fill(params_.pmeGridSaved.begin(), params_.pmeGridSaved.end(), std::complex<double>(0.0, 0.0));
    std::fill(params_.potentialGrid.begin(), params_.potentialGrid.end(), 0.0);
    
    // Allocate B-spline coefficients
    params_.bsplineCoeffs.resize(params_.splineOrder * 1000); // Assuming max 1000 atoms
}

double PGPContext::getErfcValue(double x) const {
    if (x >= Parameters::ERFC_TABLE_SIZE / Parameters::ERFC_TABLE_SCALE) {
        return 0.0;
    }
    
    int index = static_cast<int>(x * Parameters::ERFC_TABLE_SCALE);
    if (index < 0) index = 0;
    if (index >= Parameters::ERFC_TABLE_SIZE - 1) {
        return params_.erfcTable[Parameters::ERFC_TABLE_SIZE - 1];
    }
    
    // Linear interpolation
    double fraction = x * Parameters::ERFC_TABLE_SCALE - index;
    return params_.erfcTable[index] * (1.0 - fraction) + 
           params_.erfcTable[index + 1] * fraction;
}

double PGPContext::getExpValue(double x) const {
    if (x >= Parameters::EXP_TABLE_SIZE / Parameters::EXP_TABLE_SCALE) {
        return 0.0;
    }
    
    int index = static_cast<int>(x * Parameters::EXP_TABLE_SCALE);
    if (index < 0) index = 0;
    if (index >= Parameters::EXP_TABLE_SIZE - 1) {
        return params_.expTable[Parameters::EXP_TABLE_SIZE - 1];
    }
    
    // Linear interpolation
    double fraction = x * Parameters::EXP_TABLE_SCALE - index;
    return params_.expTable[index] * (1.0 - fraction) + 
           params_.expTable[index + 1] * fraction;
}

PGPContext::EnergyComponents PGPContext::computeSystemEnergy(pygcmc::model::montecarlo::MCState& state) const {
    EnergyComponents energy = {0.0, 0.0, 0.0, 0.0, 0.0};
    
    // Reset state's ewald_energy to avoid stale values
    state.ewald_energy.reset();
    
    energy.real_space = computeRealSpaceEnergy(state, false);
    energy.reciprocal = computeReciprocalEnergy(state, false);
    energy.self = computeSelfEnergy(state, false);
    energy.vdw = computeVdWEnergy(state, false);
    
    energy.total = energy.real_space + energy.reciprocal + energy.self + energy.vdw;
    
    // Store in state's ewald_energy struct for compatibility
    // Note: ewald_energy should only contain electrostatic components, not VdW
    state.ewald_energy.total = energy.real_space + energy.reciprocal + energy.self;
    state.ewald_energy.real_space = energy.real_space;
    state.ewald_energy.reciprocal = energy.reciprocal;
    state.ewald_energy.self = energy.self;
    
    return energy;
}

PGPContext::EnergyComponents PGPContext::computeMovementEnergy(
    pygcmc::model::montecarlo::MCState& state, 
    const std::vector<int>& movementResidues) const {
    
    // Suppress unused parameter warning
    (void)movementResidues;
    
    // For now, implement same as system energy
    // TODO: Implement proper movement-only calculation
    return computeSystemEnergy(state);
}

double PGPContext::computeRealSpaceEnergy(pygcmc::model::montecarlo::MCState& state, bool movement_only) const {
    // Suppress unused parameter warning
    (void)movement_only;
    
    double energy = 0.0;
    const double coulombConst = 138.935456;
    
    auto& atoms = state.atoms;
    auto& residues = state.residues;
    
    // Bounds checking
    const int maxResidues = std::min(state.activeResidueCount, static_cast<int>(residues.size()));
    const int maxAtoms = std::min(state.activeAtomCount, static_cast<int>(atoms.size()));
    
    // Double loop over residue pairs
    for (int r1 = 0; r1 < maxResidues; ++r1) {
        if (!residues[r1].active) continue;
        
        for (int r2 = r1; r2 < maxResidues; ++r2) {
            if (!residues[r2].active) continue;
            
            // Loop over atoms in each residue
            for (int i = residues[r1].atomStart; i < residues[r1].atomStart + residues[r1].atomCount; ++i) {
                if (i >= maxAtoms) break;
                
                int j_start = (r1 == r2) ? i + 1 : residues[r2].atomStart;
                for (int j = j_start; j < residues[r2].atomStart + residues[r2].atomCount; ++j) {
                    if (j >= maxAtoms) break;
                    
                    // Calculate distance with PBC
                    double dx = atoms[i].x - atoms[j].x;
                    double dy = atoms[i].y - atoms[j].y;
                    double dz = atoms[i].z - atoms[j].z;
                    
                    // Apply periodic boundary conditions
                    dx -= params_.box[0] * std::round(dx / params_.box[0]);
                    dy -= params_.box[1] * std::round(dy / params_.box[1]);
                    dz -= params_.box[2] * std::round(dz / params_.box[2]);
                    
                    double r2_val = dx*dx + dy*dy + dz*dz;
                    double r = std::sqrt(r2_val);
                    
                    if (r < params_.cutoff) {
                        // Electrostatic interaction
                        double erfc_val = getErfcValue(params_.alphaEwald * r);
                        double electrostatic = coulombConst * atoms[i].charge * atoms[j].charge * erfc_val / r;
                        
                        // Apply factor of 0.5 for double counting
                        energy += 0.5 * electrostatic;
                    }
                }
            }
        }
    }
    
    return energy;
}

double PGPContext::computeReciprocalEnergy(pygcmc::model::montecarlo::MCState& state, bool movement_only) const {
    // For PGP, reciprocal space energy comes from interpolating moveable atoms
    // in the precomputed potential grid
    
    double energy = 0.0;
    auto& atoms = state.atoms;
    auto& residues = state.residues;
    
    const int maxAtoms = std::min(state.activeAtomCount, static_cast<int>(atoms.size()));
    const int maxResidues = std::min(state.activeResidueCount, static_cast<int>(residues.size()));
    
    // Check if potential grid is properly initialized
    if (params_.potentialGrid.empty()) {
        ::pygcmc::platform::log(::pygcmc::platform::LogLevel::WARNING, "Potential grid not initialized, returning 0 reciprocal energy");
        return 0.0;
    }
    
    // Interpolate energy for moveable atoms
    for (int r = 0; r < maxResidues; ++r) {
        if (!residues[r].active) continue;
        
        // Skip fixed residues if only computing movement energy
        if (movement_only && residues[r].fixed) continue;
        
        // For standard PGP, only interpolate for non-fixed residues
        if (!movement_only && residues[r].fixed) continue;
        
        for (int i = residues[r].atomStart; 
             i < residues[r].atomStart + residues[r].atomCount; ++i) {
            if (i >= maxAtoms) break;
            
            // Skip if no charge
            if (std::abs(atoms[i].charge) < 1e-10) continue;
            
            // Calculate fractional coordinates
            double fx = atoms[i].x / params_.box[0];
            double fy = atoms[i].y / params_.box[1];
            double fz = atoms[i].z / params_.box[2];
            
            // Apply periodic boundary conditions
            fx = fx - std::floor(fx);
            fy = fy - std::floor(fy);
            fz = fz - std::floor(fz);
            
            // Scale to grid
            fx *= params_.potentialGridSize[0];
            fy *= params_.potentialGridSize[1];
            fz *= params_.potentialGridSize[2];
            
            // Simple trilinear interpolation
            int ix = static_cast<int>(fx);
            int iy = static_cast<int>(fy);
            int iz = static_cast<int>(fz);
            
            double dx = fx - ix;
            double dy = fy - iy;
            double dz = fz - iz;
            
            // Grid offset for this atom type
            int atomType = atoms[i].type;
            int gridOffset = atomType * params_.gridTotal;
            
            // Trilinear interpolation
            double potential = 0.0;
            for (int di = 0; di <= 1; ++di) {
                for (int dj = 0; dj <= 1; ++dj) {
                    for (int dk = 0; dk <= 1; ++dk) {
                        int gx = (ix + di) % params_.potentialGridSize[0];
                        int gy = (iy + dj) % params_.potentialGridSize[1];
                        int gz = (iz + dk) % params_.potentialGridSize[2];
                        
                        int gridIndex = gx * params_.gridSizeYZ + gy * params_.gridSizeZ + gz;
                        
                        double weight = (di ? dx : 1-dx) * (dj ? dy : 1-dy) * (dk ? dz : 1-dz);
                        
                        if (gridOffset + gridIndex < static_cast<int>(params_.potentialGrid.size())) {
                            potential += weight * params_.potentialGrid[gridOffset + gridIndex];
                        }
                    }
                }
            }
            
            // Energy = charge * potential
            energy += atoms[i].charge * potential;
        }
    }
    
    // The potential was halved during precomputation, so multiply by 2
    return 2.0 * energy;
}

double PGPContext::computeSelfEnergy(pygcmc::model::montecarlo::MCState& state, bool movement_only) const {
    // Suppress unused parameter warning
    (void)movement_only;
    
    double energy = 0.0;
    auto& atoms = state.atoms;
    
    const int maxAtoms = std::min(state.activeAtomCount, static_cast<int>(atoms.size()));
    
    for (int i = 0; i < maxAtoms; ++i) {
        energy += atoms[i].charge * atoms[i].charge;
    }
    
    return params_.selfEnergyCoeff * energy;
}

double PGPContext::computeVdWEnergy(pygcmc::model::montecarlo::MCState& state, bool movement_only) const {
    // Suppress unused parameter warning
    (void)movement_only;
    
    double energy = 0.0;
    auto& atoms = state.atoms;
    auto& residues = state.residues;
    auto& forcefield = state.forcefield;
    
    const int maxResidues = std::min(state.activeResidueCount, static_cast<int>(residues.size()));
    const int maxAtoms = std::min(state.activeAtomCount, static_cast<int>(atoms.size()));
    
    // Initialize residue VdW energies to zero
    for (int r = 0; r < maxResidues; ++r) {
        if (residues[r].active) {
            residues[r].energy_vdw = 0.0f;
        }
    }
    
    // Double loop over residue pairs, then atoms within residues
    for (int r1 = 0; r1 < maxResidues; ++r1) {
        if (!residues[r1].active) continue;
        
        for (int r2 = r1; r2 < maxResidues; ++r2) {
            if (!residues[r2].active) continue;
            
            double pairEnergy = 0.0;
            
            // Loop over all atom pairs between residues
            for (int a1 = residues[r1].atomStart; a1 < residues[r1].atomStart + residues[r1].atomCount; ++a1) {
                if (a1 >= maxAtoms) break;
                
                int a2_start = (r1 == r2) ? a1 + 1 : residues[r2].atomStart;
                for (int a2 = a2_start; a2 < residues[r2].atomStart + residues[r2].atomCount; ++a2) {
                    if (a2 >= maxAtoms) break;
                    
                    // Get atom types
                    int type1 = atoms[a1].type;
                    int type2 = atoms[a2].type;
                    int ljIndex = type1 * forcefield.numTotalTypes + type2;
                    
                    if (ljIndex >= static_cast<int>(forcefield.ljEps.size())) continue;
                    
                    double eps = forcefield.ljEps[ljIndex];
                    double sigma = forcefield.ljSigma[ljIndex];
                    
                    if (eps == 0.0) continue;
                    
                    // Calculate distance between atoms
                    double dx = atoms[a1].x - atoms[a2].x;
                    double dy = atoms[a1].y - atoms[a2].y;
                    double dz = atoms[a1].z - atoms[a2].z;
                    
                    // Apply PBC
                    dx -= params_.box[0] * std::round(dx / params_.box[0]);
                    dy -= params_.box[1] * std::round(dy / params_.box[1]);
                    dz -= params_.box[2] * std::round(dz / params_.box[2]);
                    
                    double r2 = dx*dx + dy*dy + dz*dz;
                    
                    if (r2 < params_.cutoff * params_.cutoff && r2 > 0.01 * 0.01) {
                        double r6 = r2 * r2 * r2;
                        double sigma6 = sigma * sigma * sigma * sigma * sigma * sigma;
                        double sigma12 = sigma6 * sigma6;
                        
                        double lj = 4.0 * eps * (sigma12 / (r6 * r6) - sigma6 / r6);
                        pairEnergy += lj;
                    }
                }
            }
            
            // Distribute energy
            if (r1 == r2) {
                // Intra-residue energy belongs entirely to this residue
                residues[r1].energy_vdw += static_cast<float>(pairEnergy);
            } else {
                // Inter-residue energy is split equally
                residues[r1].energy_vdw += static_cast<float>(pairEnergy * 0.5);
                residues[r2].energy_vdw += static_cast<float>(pairEnergy * 0.5);
            }
            energy += pairEnergy;
        }
    }
    
    return energy;
}

void PGPContext::precomputeGridPotential(pygcmc::model::montecarlo::MCState& state, int atomType) {
    // Use FFT-based approach to precompute the full electrostatic potential
    // This includes both real and reciprocal space contributions
    
    ::pygcmc::platform::log(::pygcmc::platform::LogLevel::DEBUG, "Precomputing potential grid using FFT for atom type ", atomType);
    
    // Check if we have access to global PME parameters
    if (!::pygcmc::platform::cpu::pgp_params.initialized) {
        ::pygcmc::platform::log(::pygcmc::platform::LogLevel::WARNING, "PGP global parameters not initialized, using simple real-space only");
        precomputeGridPotentialSimple(state, atomType);
        return;
    }
    
    // Save current PME grid state
    auto& pme = ::pygcmc::platform::cpu::pme_params;
    auto savedPmeGrid = pme.pmeGrid;
    auto savedMeshSize = pme.meshSize;
    
    // Temporarily adjust PME parameters to match PGP grid
    for (int i = 0; i < 3; i++) {
        pme.meshSize[i] = params_.potentialGridSize[i];
    }
    
    // Resize PME grid to match potential grid size
    int totalGridSize = params_.potentialGridSize[0] * 
                       params_.potentialGridSize[1] * 
                       params_.potentialGridSize[2];
    pme.pmeGrid.resize(totalGridSize, std::complex<double>(0.0, 0.0));
    
    // Clear the grid
    std::fill(pme.pmeGrid.begin(), pme.pmeGrid.end(), std::complex<double>(0.0, 0.0));
    
    // Spread charges from fixed atoms onto grid using PME functions
    ::pygcmc::platform::cpu::spreadChargesOntoGrid(state, true); // true = fixed only
    
    // Perform forward FFT
    ::pygcmc::platform::cpu::performFFTForward();
    
    // Apply reciprocal space factors to get potential (not energy)
    applyReciprocalSpacePotentialFactors(pme, params_);
    
    // Perform backward FFT to get potential in real space
    ::pygcmc::platform::cpu::performFFTBackward();
    
    // Apply normalization and unit conversion
    double volume = params_.box[0] * params_.box[1] * params_.box[2];
    double coulombConst = 138.935456 / params_.epsilon_r;
    double normalizationFactor = 4.0 * M_PI * coulombConst / (totalGridSize * volume);
    
    // Store the potential grid for this atom type
    int gridOffset = atomType * params_.gridTotal;
    for (int i = 0; i < totalGridSize && i < params_.gridTotal; ++i) {
        if (gridOffset + i < static_cast<int>(params_.potentialGrid.size())) {
            // Store the real part of the potential with factor of 0.5
            params_.potentialGrid[gridOffset + i] = 
                0.5 * pme.pmeGrid[i].real() * normalizationFactor;
        }
    }
    
    // Restore PME parameters
    pme.pmeGrid = savedPmeGrid;
    for (int i = 0; i < 3; i++) {
        pme.meshSize[i] = savedMeshSize[i];
    }
}

// Fallback simple implementation without FFT
void PGPContext::precomputeGridPotentialSimple(pygcmc::model::montecarlo::MCState& state, int atomType) {
    const double coulombConst = 138.935456;
    auto& atoms = state.atoms;
    auto& residues = state.residues;
    
    // Grid spacing
    double dx = params_.box[0] / params_.potentialGridSize[0];
    double dy = params_.box[1] / params_.potentialGridSize[1];
    double dz = params_.box[2] / params_.potentialGridSize[2];
    
    // Offset for this atom type in the potential grid
    int gridOffset = atomType * params_.gridTotal;
    
    // Clear the grid for this atom type
    for (int i = 0; i < params_.gridTotal; ++i) {
        if (gridOffset + i < static_cast<int>(params_.potentialGrid.size())) {
            params_.potentialGrid[gridOffset + i] = 0.0;
        }
    }
    
    // Compute potential at each grid point from fixed atoms only
    for (int gx = 0; gx < params_.potentialGridSize[0]; ++gx) {
        for (int gy = 0; gy < params_.potentialGridSize[1]; ++gy) {
            for (int gz = 0; gz < params_.potentialGridSize[2]; ++gz) {
                // Grid point position
                double x = gx * dx;
                double y = gy * dy;
                double z = gz * dz;
                
                double potential = 0.0;
                
                // Sum contribution from all fixed atoms
                const int maxResidues = std::min(state.activeResidueCount, static_cast<int>(residues.size()));
                for (int r = 0; r < maxResidues; ++r) {
                    if (!residues[r].active || !residues[r].fixed) continue;
                    
                    for (int i = residues[r].atomStart; i < residues[r].atomStart + residues[r].atomCount; ++i) {
                        if (i >= state.activeAtomCount || i >= static_cast<int>(atoms.size())) break;
                        
                        // Distance from grid point to atom
                        double dx_atom = x - atoms[i].x;
                        double dy_atom = y - atoms[i].y;
                        double dz_atom = z - atoms[i].z;
                        
                        // Apply PBC
                        dx_atom -= params_.box[0] * std::round(dx_atom / params_.box[0]);
                        dy_atom -= params_.box[1] * std::round(dy_atom / params_.box[1]);
                        dz_atom -= params_.box[2] * std::round(dz_atom / params_.box[2]);
                        
                        double r2 = dx_atom*dx_atom + dy_atom*dy_atom + dz_atom*dz_atom;
                        double r = std::sqrt(r2);
                        
                        if (r > 0.01 && r < params_.potential_cutoff) {
                            // Electrostatic potential
                            double erfc_val = getErfcValue(params_.alphaEwald * r);
                            potential += coulombConst * atoms[i].charge * erfc_val / r;
                        }
                    }
                }
                
                // Store in grid with factor of 0.5
                int gridIndex = gx * params_.gridSizeYZ + gy * params_.gridSizeZ + gz;
                if (gridOffset + gridIndex < static_cast<int>(params_.potentialGrid.size())) {
                    params_.potentialGrid[gridOffset + gridIndex] = 0.5 * potential;
                }
            }
        }
    }
}

// Apply reciprocal space factors to convert from charge density to potential
void PGPContext::applyReciprocalSpacePotentialFactors(::pygcmc::platform::cpu::PMEParams& pme, const PGPContext::Parameters& params) {
    int nx = pme.meshSize[0];
    int ny = pme.meshSize[1];
    int nz = pme.meshSize[2];
    
    double factor = M_PI * M_PI / (params.alphaEwald * params.alphaEwald);
    
    int maxkx = (nx + 1) / 2;
    int maxky = (ny + 1) / 2;
    int maxkz = (nz + 1) / 2;
    
    // Calculate reciprocal lattice vectors
    double recipBoxVectors[3][3] = {{0}};
    recipBoxVectors[0][0] = 2.0 * M_PI / params.box[0];
    recipBoxVectors[1][1] = 2.0 * M_PI / params.box[1];
    recipBoxVectors[2][2] = 2.0 * M_PI / params.box[2];
    
    for (int kx = 0; kx < nx; kx++) {
        double mx = (kx < maxkx) ? kx : (kx - nx);
        double mhx = mx * recipBoxVectors[0][0];
        
        for (int ky = 0; ky < ny; ky++) {
            double my = (ky < maxky) ? ky : (ky - ny);
            double mhy = my * recipBoxVectors[1][1];
            
            for (int kz = 0; kz < nz; kz++) {
                if (kx == 0 && ky == 0 && kz == 0) continue;
                
                double mz = (kz < maxkz) ? kz : (kz - nz);
                double mhz = mz * recipBoxVectors[2][2];
                
                double k2 = mhx * mhx + mhy * mhy + mhz * mhz;
                
                // Get B-spline moduli
                double bx = pme.bsplineModuli[0][kx];
                double by = pme.bsplineModuli[1][ky];
                double bz = pme.bsplineModuli[2][kz];
                double denom = k2 * bx * by * bz;
                
                if (denom < 1e-10) denom = 1e-10;
                
                // Apply factor for potential: exp(-k²/(4α²)) / (k² * B)
                double potentialFactor = std::exp(-k2 * factor / 4.0) / denom;
                
                int index = kx * ny * nz + ky * nz + kz;
                pme.pmeGrid[index] *= potentialFactor;
            }
        }
    }
}

} // namespace pgp
} // namespace energy
} // namespace cpu
} // namespace platform