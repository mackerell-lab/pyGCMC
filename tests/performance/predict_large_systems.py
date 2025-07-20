#!/usr/bin/env python
"""Predict performance for large water systems based on measured scaling"""

import numpy as np

# Measured scaling law: time = 4.926e-05 * N^1.96
a = 4.926e-05
b = 1.96

print("Performance Predictions for Large Water Systems")
print("Based on measured scaling: time = 4.926e-05 * N^1.96")
print("=" * 70)

# System sizes to predict
sizes = [
    (7, 343),
    (8, 512), 
    (10, 1000),
    (15, 3375),
    (20, 8000),
    (25, 15625),
    (30, 27000),
    (32, 32768),  # ~10 nm box
]

print(f"{'n/dim':>5} | {'N waters':>8} | {'Box (nm)':>8} | {'ms/step':>10} | {'steps/s':>10} | {'20k steps':>12}")
print("-" * 70)

for n_dim, n_waters in sizes:
    box_size = n_dim * 0.31
    time_ms = a * (n_waters ** b)
    steps_per_sec = 1000 / time_ms
    time_20k = time_ms * 20000 / 1000  # seconds
    
    # Format time nicely
    if time_20k < 60:
        time_str = f"{time_20k:.1f} s"
    elif time_20k < 3600:
        time_str = f"{time_20k/60:.1f} min"
    else:
        time_str = f"{time_20k/3600:.1f} hr"
    
    print(f"{n_dim:5d} | {n_waters:8d} | {box_size:8.1f} | {time_ms:10.2f} | {steps_per_sec:10.1f} | {time_str:>12}")

print("\n" + "=" * 70)
print("PRACTICAL SYSTEM SIZES")
print("=" * 70)

# Find practical limits
target_times = [1, 10, 60, 600]  # seconds for 20k steps
print("\nMaximum system size for different time targets (20,000 steps):")

for target_sec in target_times:
    target_ms = target_sec * 1000 / 20000
    # Solve: target_ms = a * N^b for N
    max_n = (target_ms / a) ** (1/b)
    n_dim = int((max_n ** (1/3)) + 0.5)
    actual_n = n_dim ** 3
    actual_time = a * (actual_n ** b) * 20000 / 1000
    
    if target_sec < 60:
        target_str = f"{target_sec} sec"
    else:
        target_str = f"{target_sec/60:.0f} min"
    
    print(f"  {target_str:>6}: up to {actual_n:5d} waters ({n_dim}×{n_dim}×{n_dim}, {n_dim*0.31:.1f} nm box) → {actual_time:.1f} sec")

print("\n" + "=" * 70)
print("DRUDE MODEL PREDICTIONS")
print("=" * 70)
print("\nAssuming Drude is 10-50x slower due to SCF iterations:")

for factor in [10, 20, 50]:
    print(f"\nIf Drude is {factor}x slower:")
    for n_dim, n_waters in [(5, 125), (6, 216), (7, 343), (8, 512)]:
        time_ms = a * (n_waters ** b) * factor
        time_20k = time_ms * 20000 / 1000
        if time_20k < 3600:
            print(f"  {n_waters:4d} waters: {time_20k/60:5.1f} minutes")
        else:
            print(f"  {n_waters:4d} waters: {time_20k/3600:5.1f} hours")