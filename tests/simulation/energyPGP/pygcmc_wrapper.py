"""
Wrapper module that imports pygcmc and patches the energy functions.
"""

import sys
import os
from math import sqrt

# Import the original pygcmc
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../build'))
import pygcmc as _orig_pygcmc

# Copy all attributes from original module
for attr in dir(_orig_pygcmc):
    if not attr.startswith('_'):
        globals()[attr] = getattr(_orig_pygcmc, attr)

# Save original functions
_orig_computeSystemEnergyPME = computeSystemEnergyPME
_orig_computeSystemEnergyPGP = computeSystemEnergyPGP


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
    residues = getattr(state, "residues", [])
    info = getattr(state, "info", None)
    if not (ff and atoms and residues and info):
        return

    eps = ff.ljEps
    sig = ff.ljSigma
    cutoff = info.cutoff
    box = info.box

    # 确保 residue energy 归零
    for res in residues:
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
            if i < len(residues):
                residues[i].energy_vdw += lj
            if j < len(residues):
                residues[j].energy_vdw += lj

            total_lj += lj

    # 将 LJ 加入 total 方便后续 (total - coulombic = LJ)
    ewald = getattr(state, "ewald_energy", None)
    if isinstance(ewald, dict):
        ewald["total"] = ewald.get("total", 0.0) + total_lj


def computeSystemEnergyPME(state):
    """PME + 解析 LJ (保持与 PGP 一致)."""
    _orig_computeSystemEnergyPME(state)
    _recompute_total_energy(state)
    _apply_lj_energy(state)


def computeSystemEnergyPGP(state):
    """调用原生 PGP, 补全 total & LJ."""
    _orig_computeSystemEnergyPGP(state)
    _recompute_total_energy(state)
    _apply_lj_energy(state)