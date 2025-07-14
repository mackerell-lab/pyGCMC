#include "PGPCoreIndependent.hpp"
#include "platform/platform.hpp"
#include "../common/EnergyConstants.hpp"
#include <cmath>
#include <algorithm>

namespace pygcmc {
namespace platform {
namespace cpu {

// Global independent PGP parameters instance
PGPParamsIndependent pgp_params_independent;

// Mathematical constants
const double TWO_PI = 2.0 * M_PI;
const double SQRT_PI = sqrt(M_PI);

// Forward declarations for helper functions
void interpolateMoleculeEnergyIndependent(model::MCState& state, double& energy, bool movement_only = false);
void computeRealSpacePGPIndependent(model::MCState& state, bool movement_only);
double computeSelfEnergyPGPIndependent(model::MCState& state, bool movement_only);
void precomputeGridPotentialIndependent(model::MCState& state, bool fixed_only);

void PGPParamsIndependent::initializeOwnTables() {
    // Initialize erfc and Ewald scaling tables independently
    const int tableSize = 2000;
    const double maxR = cutoff * 1.2;  // Extend slightly beyond cutoff
    
    erfcTable.resize(tableSize);
    ewaldScaleTable.resize(tableSize);
    
    ewaldDX = maxR / (tableSize - 1);
    ewaldDXInv = 1.0 / ewaldDX;
    erfcDXInv = ewaldDXInv;  // Same spacing for both tables
    
    // Fill tables
    for (int i = 0; i < tableSize; ++i) {
        double r = i * ewaldDX;
        double alphaR = alpha * r;
        
        // erfc table
        erfcTable[i] = std::erfc(alphaR);
        
        // Ewald scaling table: exp(-alpha^2 * r^2)
        ewaldScaleTable[i] = std::exp(-alphaR * alphaR);
    }
    
    platform::log(LogLevel::DEBUG, "PGP: Initialized independent lookup tables");
}

void PGPParamsIndependent::initializeOwnBsplines() {
    // Initialize B-spline moduli independently
    for (int dim = 0; dim < 3; ++dim) {
        int meshPoints = meshSize[dim];
        bsplineModuli[dim].resize(meshPoints);
        
        // Calculate B-spline moduli for each k-vector
        for (int i = 0; i < meshPoints; ++i) {
            double k = 2.0 * M_PI * i / meshPoints;
            if (i > meshPoints / 2) {
                k -= 2.0 * M_PI;
            }
            
            // B-spline modulus calculation
            double sum = 0.0;
            for (int j = 0; j < splineOrder; ++j) {
                double arg = k * j / splineOrder;
                sum += std::cos(arg);
            }
            
            double modulus = std::pow(sum / splineOrder, splineOrder);
            bsplineModuli[dim][i] = (modulus > 1e-10) ? 1.0 / modulus : 0.0;
        }
    }
    
    platform::log(LogLevel::DEBUG, "PGP: Initialized independent B-spline moduli");
}

void PGPParamsIndependent::initializeOwnGrids() {
    // Initialize FFT grids independently
    int totalGridPoints = meshSize[0] * meshSize[1] * meshSize[2];
    
    // Allocate own grids
    pmeGrid.resize(totalGridPoints);
    pmeCharge.resize(totalGridPoints);
    
    // Clear grids
    std::fill(pmeGrid.begin(), pmeGrid.end(), std::complex<double>(0.0, 0.0));
    std::fill(pmeCharge.begin(), pmeCharge.end(), 0.0);
    
    platform::log(LogLevel::DEBUG, "PGP: Allocated independent FFT grids: ", 
                  meshSize[0], "x", meshSize[1], "x", meshSize[2]);
}

void PGPParamsIndependent::initializePotentialGrid() {
    // Initialize potential grid for PGP
    int totalPotentialPoints = potential_grid_size[0] * 
                              potential_grid_size[1] * 
                              potential_grid_size[2];
    
    potentialGrid.resize(totalPotentialPoints);
    std::fill(potentialGrid.begin(), potentialGrid.end(), std::complex<double>(0.0, 0.0));
    
    // Calculate grid spacing
    grid_spacing = potential_cutoff / std::min({potential_grid_size[0], 
                                                potential_grid_size[1], 
                                                potential_grid_size[2]});
    
    platform::log(LogLevel::DEBUG, "PGP: Initialized potential grid: ",
                  potential_grid_size[0], "x", potential_grid_size[1], "x", 
                  potential_grid_size[2], " spacing=", grid_spacing);
}

void setPGPParametersIndependent(double alpha, const int meshSize[3], 
                                double potential_cutoff, 
                                const int potentialGridSize[3], 
                                int splineOrder, double tolerance) {
    // Set parameters without touching PME
    pgp_params_independent.alpha = alpha;
    pgp_params_independent.tolerance = tolerance;
    pgp_params_independent.epsilon_r = 1.0;  // Default
    pgp_params_independent.splineOrder = splineOrder;
    pgp_params_independent.potential_cutoff = potential_cutoff;
    
    // Copy array parameters
    for (int i = 0; i < 3; i++) {
        pgp_params_independent.meshSize[i] = meshSize[i];
        pgp_params_independent.potential_grid_size[i] = potentialGridSize[i];
    }
    
    platform::log(LogLevel::INFO, "PGP: Set independent parameters");
}

void initializePGPParametersIndependent(double cutoff, const double box[3], 
                                       double alpha, const int meshSize[3], 
                                       double potentialCutoff,
                                       const int potentialGridSize[3],
                                       int splineOrder, double tolerance) {
    // Set basic parameters
    pgp_params_independent.cutoff = cutoff;
    for (int i = 0; i < 3; i++) {
        pgp_params_independent.box[i] = box[i];
    }
    
    // Set PGP-specific parameters
    setPGPParametersIndependent(alpha, meshSize, potentialCutoff, 
                               potentialGridSize, splineOrder, tolerance);
    
    // Initialize all data structures independently
    pgp_params_independent.initializeIndependent();
    
    platform::log(LogLevel::INFO, "PGP: Initialized with independent data structures");
}

void computeSystemEnergyPGPIndependent(model::MCState& state) {
    // Implementation using pgp_params_independent instead of shared params
    if (!pgp_params_independent.initialized) {
        throw std::runtime_error("PGP independent parameters not initialized");
    }
    
    platform::log(LogLevel::DEBUG, "Computing total system energy using independent PGP method");
    
    // 1. Calculate grid potential interpolation part
    double grid_energy = 0.0;
    interpolateMoleculeEnergyIndependent(state, grid_energy, false);
    
    // 2. Calculate real space part
    computeRealSpacePGPIndependent(state, false);
    
    // 3. Calculate self energy correction
    state.ewald_energy.self = computeSelfEnergyPGPIndependent(state, false);
    
    // 4. Calculate LJ interactions would be done separately
    // Since we're focusing on electrostatics, set to 0
    double vdw_total = 0.0;
    
    // Multiply real space energy by COULOMB constant
    state.ewald_energy.real_space *= COULOMB;
    
    // Calculate total energy
    state.ewald_energy.reciprocal = grid_energy;
    state.ewald_energy.total = grid_energy + state.ewald_energy.real_space + 
                             state.ewald_energy.self + vdw_total;
    
    platform::log(LogLevel::DEBUG, "PGP independent system energy components: ",
                  "Grid=", grid_energy, 
                  ", Real=", state.ewald_energy.real_space,
                  ", Self=", state.ewald_energy.self,
                  ", Total=", state.ewald_energy.total);
}

void computeMovementEnergyPGPIndependent(model::MCState& state) {
    // Implementation using pgp_params_independent
    if (!pgp_params_independent.initialized) {
        throw std::runtime_error("PGP independent parameters not initialized");
    }
    
    platform::log(LogLevel::DEBUG, "Computing movement residue energy using independent PGP method");
    
    // 1. Calculate grid potential interpolation part
    double grid_energy = 0.0;
    interpolateMoleculeEnergyIndependent(state, grid_energy, true);
    
    // 2. Calculate real space part (only for moving residues)
    computeRealSpacePGPIndependent(state, true);
    
    // 3. Calculate self energy correction (only for moving residues)
    state.ewald_energy.self = computeSelfEnergyPGPIndependent(state, true);
    
    // 4. Calculate LJ interactions would be done separately
    double vdw_total = 0.0;
    
    // Multiply real space energy by COULOMB constant
    state.ewald_energy.real_space *= COULOMB;
    
    // Calculate total energy
    state.ewald_energy.reciprocal = grid_energy;
    state.ewald_energy.total = grid_energy + state.ewald_energy.real_space + 
                             state.ewald_energy.self + vdw_total;
    
    platform::log(LogLevel::DEBUG, "PGP independent movement energy components: ",
                  "Grid=", grid_energy,
                  ", Real=", state.ewald_energy.real_space,
                  ", Self=", state.ewald_energy.self,
                  ", Total=", state.ewald_energy.total);
}

// Helper function implementations

void interpolateMoleculeEnergyIndependent(model::MCState& state, double& energy, bool movement_only) {
    // Reset energy
    energy = 0.0;
    double raw_energy = 0.0;
    
    // Use actual vector sizes for safety
    const int maxResidues = std::min(state.activeResidueCount, static_cast<int>(state.residues.size()));
    const int maxAtoms = std::min(state.activeAtomCount, static_cast<int>(state.atoms.size()));
    
    // Lambda function to process a residue
    auto process_residue = [&](int res_idx) {
        if (res_idx < 0 || res_idx >= maxResidues) return;
        const auto& residue = state.residues[res_idx];
        
        if (!residue.active || residue.fixed) return;
        
        // Process atoms in this residue
        const int atom_start = residue.atomStart;
        const int atom_end = std::min(atom_start + residue.atomCount, maxAtoms);
        
        for (int j = atom_start; j < atom_end; j++) {
            if (j < 0 || j >= maxAtoms) continue;
            const auto& atom = state.atoms[j];
            
            // Only process charged atoms
            if (std::abs(atom.charge) < 1e-6) continue;
            
            // Calculate grid position
            double pos[3] = {atom.x, atom.y, atom.z};
            
            // Calculate fractional coordinates
            double fractional[3];
            for (int d = 0; d < 3; d++) {
                fractional[d] = pos[d] / pgp_params_independent.box[d];
                fractional[d] -= floor(fractional[d]);  // Ensure in [0,1) range
                fractional[d] *= pgp_params_independent.potential_grid_size[d];
            }
            
            // Calculate grid index and fractional part
            int gridIndices[3];
            double gridFractions[3];
            for (int d = 0; d < 3; d++) {
                gridFractions[d] = fractional[d] - floor(fractional[d]);
                gridIndices[d] = static_cast<int>(floor(fractional[d]));
                if (gridIndices[d] < 0) 
                    gridIndices[d] += pgp_params_independent.potential_grid_size[d];
            }
            
            // Simple trilinear interpolation for now
            int nx = pgp_params_independent.potential_grid_size[0];
            int ny = pgp_params_independent.potential_grid_size[1];
            int nz = pgp_params_independent.potential_grid_size[2];
            
            double potential = 0.0;
            
            // 8-point interpolation
            for (int ix = 0; ix < 2; ix++) {
                int xindex = (gridIndices[0] + ix) % nx;
                double wx = (ix == 0) ? (1.0 - gridFractions[0]) : gridFractions[0];
                
                for (int iy = 0; iy < 2; iy++) {
                    int yindex = (gridIndices[1] + iy) % ny;
                    double wy = (iy == 0) ? (1.0 - gridFractions[1]) : gridFractions[1];
                    
                    for (int iz = 0; iz < 2; iz++) {
                        int zindex = (gridIndices[2] + iz) % nz;
                        double wz = (iz == 0) ? (1.0 - gridFractions[2]) : gridFractions[2];
                        
                        int index = xindex * ny * nz + yindex * nz + zindex;
                        double grid_value = pgp_params_independent.potentialGrid[index].real();
                        potential += grid_value * wx * wy * wz;
                    }
                }
            }
            
            // Accumulate energy
            raw_energy += potential * atom.charge;
        }
    };
    
    if (movement_only) {
        // Only process movement residues
        for (const auto& movementInfo : state.movementResidues) {
            const int start = movementInfo.startIndex;
            const int end = std::min(start + movementInfo.activeCount, maxResidues);
            
            for (int i = start; i < end; i++) {
                process_residue(i);
            }
        }
    } else {
        // Process all non-fixed active residues
        for (int i = 0; i < maxResidues; i++) {
            process_residue(i);
        }
    }
    
    // Factor of 2 for PGP convention
    energy = 2.0 * raw_energy;
}

void computeRealSpacePGPIndependent(model::MCState& state, bool movement_only) {
    const auto& box = state.info.box;
    auto& atoms = state.atoms;
    auto& residues = state.residues;
    const double cutoff2 = pgp_params_independent.cutoff * pgp_params_independent.cutoff;
    
    // Reset real space energy
    state.ewald_energy.real_space = 0.0;
    
    // Use actual vector sizes for safety
    const int maxResidues = std::min(state.activeResidueCount, static_cast<int>(residues.size()));
    const int maxAtoms = std::min(state.activeAtomCount, static_cast<int>(atoms.size()));
    
    // Loop over all residue pairs
    for(int r1 = 0; r1 < maxResidues; r1++) {
        if(!residues[r1].active) continue;
        
        bool r1_in_movement = false;
        if(movement_only) {
            for(const auto& movementInfo : state.movementResidues) {
                if(r1 >= movementInfo.startIndex && 
                   r1 < movementInfo.startIndex + movementInfo.activeCount) {
                    r1_in_movement = true;
                    break;
                }
            }
        }
        
        for(int r2 = 0; r2 < maxResidues; r2++) {
            if(r1 == r2) continue;  // Skip self-interaction
            if(!residues[r2].active) continue;
            
            // For movement_only mode, at least one residue must be moving
            if(movement_only) {
                bool r2_in_movement = false;
                for(const auto& movementInfo : state.movementResidues) {
                    if(r2 >= movementInfo.startIndex && 
                       r2 < movementInfo.startIndex + movementInfo.activeCount) {
                        r2_in_movement = true;
                        break;
                    }
                }
                if(!r1_in_movement && !r2_in_movement) continue;
            }
            
            // Skip if both are fixed (already in grid)
            if(residues[r1].fixed && residues[r2].fixed) continue;
            
            // To avoid double counting, only calculate once per pair
            if(r2 < r1) continue;
            
            // Loop over atoms with bounds checking
            const int i_start = residues[r1].atomStart;
            const int i_end = std::min(i_start + residues[r1].atomCount, maxAtoms);
            
            const int j_start = residues[r2].atomStart;
            const int j_end = std::min(j_start + residues[r2].atomCount, maxAtoms);
            
            for(int i = i_start; i < i_end; i++) {
                if(i < 0 || i >= maxAtoms) continue;
                
                for(int j = j_start; j < j_end; j++) {
                    if(j < 0 || j >= maxAtoms) continue;
                    
                    float dx = atoms[i].x - atoms[j].x;
                    float dy = atoms[i].y - atoms[j].y;
                    float dz = atoms[i].z - atoms[j].z;
                    
                    // Apply PBC
                    dx -= box[0] * round(dx / box[0]);
                    dy -= box[1] * round(dy / box[1]);
                    dz -= box[2] * round(dz / box[2]);
                    
                    float r2 = dx*dx + dy*dy + dz*dz;
                    
                    if(r2 > cutoff2) continue;
                    
                    float r = sqrt(r2);
                    float qi = atoms[i].charge;
                    float qj = atoms[j].charge;
                    
                    if(std::abs(qi) < 1e-6 || std::abs(qj) < 1e-6) continue;
                    
                    // Use independent erfc approximation
                    double alphaR = pgp_params_independent.alpha * r;
                    double term = std::erfc(alphaR);
                    double pair_energy = qi * qj * term / r;
                    
                    state.ewald_energy.real_space += pair_energy;
                }
            }
        }
    }
}

double computeSelfEnergyPGPIndependent(model::MCState& state, bool movement_only) {
    double self_energy = 0.0;
    const double factor = -pgp_params_independent.alpha / SQRT_PI;
    
    // Use actual vector sizes for safety
    const int maxResidues = std::min(state.activeResidueCount, static_cast<int>(state.residues.size()));
    const int maxAtoms = std::min(state.activeAtomCount, static_cast<int>(state.atoms.size()));
    
    if (movement_only) {
        // Only sum charges for movement residues
        for (const auto& movementInfo : state.movementResidues) {
            const int start = movementInfo.startIndex;
            const int end = std::min(start + movementInfo.activeCount, maxResidues);
            
            for (int i = start; i < end; i++) {
                if (i < 0 || i >= maxResidues) continue;
                if (!state.residues[i].active) continue;
                
                const int atom_start = state.residues[i].atomStart;
                const int atom_end = std::min(atom_start + state.residues[i].atomCount, maxAtoms);
                
                for (int j = atom_start; j < atom_end; j++) {
                    if (j < 0 || j >= maxAtoms) continue;
                    double charge = state.atoms[j].charge;
                    self_energy += charge * charge;
                }
            }
        }
    } else {
        // Sum all active residue charges
        for (int i = 0; i < maxResidues; i++) {
            if (!state.residues[i].active) continue;
            
            const int atom_start = state.residues[i].atomStart;
            const int atom_end = std::min(atom_start + state.residues[i].atomCount, maxAtoms);
            
            for (int j = atom_start; j < atom_end; j++) {
                if (j < 0 || j >= maxAtoms) continue;
                double charge = state.atoms[j].charge;
                self_energy += charge * charge;
            }
        }
    }
    
    return self_energy * factor * COULOMB;
}

} // namespace cpu
} // namespace platform
} // namespace pygcmc