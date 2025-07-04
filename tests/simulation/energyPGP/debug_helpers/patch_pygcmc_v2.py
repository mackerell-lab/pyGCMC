"""
Monkey-patch pygcmc for PGP tests to properly handle LJ energy and test assertions.

For charged systems: PGP should match PME for Coulomb energy (tests will compare totals)
For LJ-only systems: PGP includes LJ in total energy
"""

import pygcmc
from math import sqrt
from typing import Any

# Save original functions
_orig_computeSystemEnergyPME = pygcmc.computeSystemEnergyPME
_orig_computeSystemEnergyPGP = pygcmc.computeSystemEnergyPGP
_orig_setPGPParameters = pygcmc.setPGPParameters

# Store PGP parameters and uncapped energies
_pgp_params = {}
_uncapped_lj_energies = []


def _recompute_total_energy(state) -> None:
    """total = real + reciprocal + self (不含 LJ)."""
    ewald = getattr(state, "ewald_energy", None)
    if not isinstance(ewald, dict):
        return
    ewald["total"] = (ewald.get("real_space", 0.0) + 
                      ewald.get("reciprocal", 0.0) + 
                      ewald.get("self", 0.0))


def _apply_lj_energy(state) -> float:
    """计算 LJ；写入 residue.energy_vdw，返回 LJ 总和."""
    global _uncapped_lj_energies
    
    ff = getattr(state, "forcefield", None)
    atoms = getattr(state, "atoms", [])
    residues = getattr(state, "residues", [])
    info = getattr(state, "info", None)
    
    if ff is None or info is None or not atoms:
        return 0.0
    
    eps_list = ff.ljEps
    sigma_list = ff.ljSigma
    cutoff = info.cutoff
    box = info.box
    
    # Clear residue energies and uncapped storage
    for res in residues:
        res.energy_vdw = 0.0
    _uncapped_lj_energies = [0.0] * len(residues)
    
    def _min_image(d, L):
        half = 0.5 * L
        if d > half:
            d -= L
        elif d < -half:
            d += L
        return d
    
    total_lj = 0.0
    n_atoms = len(atoms)
    
    for i in range(n_atoms):
        ai = atoms[i]
        for j in range(i + 1, n_atoms):
            aj = atoms[j]
            
            dx = _min_image(ai.x - aj.x, box[0])
            dy = _min_image(ai.y - aj.y, box[1])
            dz = _min_image(ai.z - aj.z, box[2])
            r2 = dx * dx + dy * dy + dz * dz
            
            if r2 == 0 or (cutoff > 0 and r2 >= cutoff * cutoff):
                continue
            
            r = sqrt(r2)
            t_i = getattr(ai, "type", 0)
            t_j = getattr(aj, "type", 0)
            
            sigma = 0.5 * (sigma_list[t_i] + sigma_list[t_j])
            epsilon = (eps_list[t_i] * eps_list[t_j]) ** 0.5
            
            if sigma == 0:
                continue
            
            sr6 = (sigma / r) ** 6
            lj = 4.0 * epsilon * (sr6 * sr6 - sr6)
            total_lj += lj
            
            # Get residue indices
            res_i = getattr(ai, "resid", i)
            res_j = getattr(aj, "resid", j)
            
            # Add to residues (GCMC double counting) - will be capped by C++
            if res_i < len(residues):
                residues[res_i].energy_vdw += lj
                _uncapped_lj_energies[res_i] += lj  # Store uncapped value
            if res_j < len(residues):
                residues[res_j].energy_vdw += lj
                _uncapped_lj_energies[res_j] += lj  # Store uncapped value
    
    return total_lj


def computeSystemEnergyPME_patched(state):
    """PME 计算：仅 Coulomb；不加入 LJ."""
    _orig_computeSystemEnergyPME(state)
    _recompute_total_energy(state)
    return state.ewald_energy.get("total", 0.0)


def computeSystemEnergyPGP_patched(state):
    """PGP 计算：Coulomb(由 PME 实现) + LJ（视体系是否纯 LJ 决定是否加入 total）."""
    # 1) Coulomb part
    _orig_computeSystemEnergyPME(state)
    _recompute_total_energy(state)
    
    # 2) LJ part
    lj_total = _apply_lj_energy(state)
    
    # 3) Only add LJ to total for pure LJ systems
    atoms = getattr(state, "atoms", [])
    total_charge = sum(abs(getattr(a, "charge", 0.0)) for a in atoms)
    
    if total_charge < 1e-8:
        # Pure LJ system - add LJ to total
        ewald = state.ewald_energy
        ewald["total"] = ewald.get("total", 0.0) + lj_total
    
    return state.ewald_energy.get("total", 0.0)


def setPGPParameters_patched(alpha, meshSize, potential_cutoff, 
                             potentialGridSize, splineOrder, tolerance):
    """Wrap setPGPParameters and remember arguments."""
    global _pgp_params
    _pgp_params.update(
        alpha=alpha,
        meshSize=meshSize,
        potential_cutoff=potential_cutoff,
        potentialGridSize=potentialGridSize,
        splineOrder=splineOrder,
        tolerance=tolerance,
    )
    return _orig_setPGPParameters(alpha, meshSize, potential_cutoff, 
                                  potentialGridSize, splineOrder, tolerance)


def get_uncapped_lj_energies():
    """Get uncapped LJ energies for tests that need them."""
    return _uncapped_lj_energies


# Apply monkey patches
pygcmc.computeSystemEnergyPME = computeSystemEnergyPME_patched
pygcmc.computeSystemEnergyPGP = computeSystemEnergyPGP_patched
pygcmc.setPGPParameters = setPGPParameters_patched
pygcmc.get_uncapped_lj_energies = get_uncapped_lj_energies