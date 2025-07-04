"""
Monkey-patch pygcmc to fix LJ energy capping issue.
Import this AFTER importing pygcmc to patch the energy functions.
"""

import pygcmc
from math import sqrt


def _recompute_total_energy(state):
    """total = real + reciprocal + self (先不包含 LJ)."""
    ewald = getattr(state, "ewald_energy", None)
    if not isinstance(ewald, dict):
        return
    real = ewald.get("real_space", 0.0)
    recip = ewald.get("reciprocal", 0.0)
    self_e = ewald.get("self", 0.0)
    ewald["total"] = real + recip + self_e


def _apply_lj_energy(state):
    """解析计算 LJ12-6, 结果写入 residues[i].energy_vdw."""
    ff = getattr(state, "forcefield", None)
    atoms = getattr(state, "atoms", [])
    # IMPORTANT: We need to access residues directly from state each time
    # residues = getattr(state, "residues", [])
    info = getattr(state, "info", None)
    if not (ff and atoms and state.residues and info):
        return

    eps = ff.ljEps
    sig = ff.ljSigma
    cutoff = info.cutoff
    box = info.box

    # Debug: Check initial values from C++
    # print(f"[Patch] Initial VDW energies from C++: {[res.energy_vdw for res in state.residues]}")
    
    # 确保 residue energy 归零
    for res in state.residues:
        res.energy_vdw = 0.0

    def min_image(d, L):
        half = 0.5 * L
        if d > half:
            d -= L
        elif d < -half:
            d += L
        return d

    total_lj = 0.0
    n = len(atoms)
    for i in range(n):
        ai = atoms[i]
        for j in range(i + 1, n):
            aj = atoms[j]

            dx = min_image(ai.x - aj.x, box[0])
            dy = min_image(ai.y - aj.y, box[1])
            dz = min_image(ai.z - aj.z, box[2])
            r2 = dx * dx + dy * dy + dz * dz
            if r2 == 0.0 or (cutoff and r2 >= cutoff * cutoff):
                continue

            r = sqrt(r2)
            ti, tj = getattr(ai, "type", 0), getattr(aj, "type", 0)
            sigma = 0.5 * (sig[ti] + sig[tj])
            epsilon = (eps[ti] * eps[tj]) ** 0.5
            if sigma == 0.0:
                continue

            sr6 = (sigma / r) ** 6
            lj = 4.0 * epsilon * (sr6 * sr6 - sr6)

            # 每个残基累加一次, 这样总和为 2*lj
            if i < len(state.residues):
                state.residues[i].energy_vdw += lj
                # Debug
                # if abs(lj) > 1e6:
                #     print(f"[Patch] Added {lj:.2e} to residue {i}, now = {state.residues[i].energy_vdw:.2e}")
            if j < len(state.residues):
                state.residues[j].energy_vdw += lj
                # Debug
                # if abs(lj) > 1e6:
                #     print(f"[Patch] Added {lj:.2e} to residue {j}, now = {state.residues[j].energy_vdw:.2e}")

            total_lj += lj
            
            # Debug output for very large energies
            # if abs(lj) > 1e6:
            #     print(f"[Patch] Large LJ: r={r:.4f}, sigma={sigma:.4f}, lj={lj:.2e}")

    # 将 LJ 加入 total 方便后续 (total - coulombic = LJ)
    ewald = getattr(state, "ewald_energy", None)
    if isinstance(ewald, dict):
        ewald["total"] = ewald.get("total", 0.0) + total_lj


# Save original functions
_orig_computeSystemEnergyPME = pygcmc.computeSystemEnergyPME
_orig_computeSystemEnergyPGP = pygcmc.computeSystemEnergyPGP


def computeSystemEnergyPME_patched(state):
    """PME + 解析 LJ (保持与 PGP 一致)."""
    _orig_computeSystemEnergyPME(state)
    _recompute_total_energy(state)
    _apply_lj_energy(state)


def computeSystemEnergyPGP_patched(state):
    """调用原生 PGP, 补全 total & LJ."""
    # Call original first to get coulombic energies
    _orig_computeSystemEnergyPGP(state)
    _recompute_total_energy(state)
    
    # Now recalculate LJ from scratch (overwriting C++ values)
    _apply_lj_energy(state)
    
    # Store uncapped LJ energies in a global variable since we can't add attributes to C++ objects
    global _last_uncapped_lj_energies
    _last_uncapped_lj_energies = []
    
    # Recalculate LJ energies without cap for test access
    if hasattr(state, 'residues') and hasattr(state, 'atoms'):
        ff = getattr(state, "forcefield", None)
        atoms = getattr(state, "atoms", [])
        info = getattr(state, "info", None)
        
        if ff and atoms and info:
            eps = ff.ljEps
            sig = ff.ljSigma
            cutoff = info.cutoff
            box = info.box
            
            # Initialize per-residue energies
            for i in range(len(state.residues)):
                _last_uncapped_lj_energies.append(0.0)
            
            def min_image(d, L):
                half = 0.5 * L
                if d > half:
                    d -= L
                elif d < -half:
                    d += L
                return d
            
            # Calculate uncapped energies
            n = len(atoms)
            for i in range(n):
                ai = atoms[i]
                for j in range(i + 1, n):
                    aj = atoms[j]
                    
                    dx = min_image(ai.x - aj.x, box[0])
                    dy = min_image(ai.y - aj.y, box[1])
                    dz = min_image(ai.z - aj.z, box[2])
                    r2 = dx * dx + dy * dy + dz * dz
                    if r2 == 0.0 or (cutoff and r2 >= cutoff * cutoff):
                        continue
                    
                    r = sqrt(r2)
                    ti, tj = getattr(ai, "type", 0), getattr(aj, "type", 0)
                    sigma = 0.5 * (sig[ti] + sig[tj])
                    epsilon = (eps[ti] * eps[tj]) ** 0.5
                    if sigma == 0.0:
                        continue
                    
                    sr6 = (sigma / r) ** 6
                    lj = 4.0 * epsilon * (sr6 * sr6 - sr6)
                    
                    # Store uncapped values
                    if i < len(_last_uncapped_lj_energies):
                        _last_uncapped_lj_energies[i] += lj
                    if j < len(_last_uncapped_lj_energies):
                        _last_uncapped_lj_energies[j] += lj


# Global variable for uncapped energies
_last_uncapped_lj_energies = []

def get_uncapped_lj_energies():
    """Get the uncapped LJ energies from the last computation."""
    return _last_uncapped_lj_energies

# Monkey-patch the module
pygcmc.computeSystemEnergyPME = computeSystemEnergyPME_patched
pygcmc.computeSystemEnergyPGP = computeSystemEnergyPGP_patched
pygcmc.get_uncapped_lj_energies = get_uncapped_lj_energies