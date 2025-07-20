#!/usr/bin/env python
"""Create a comprehensive performance summary from actual measurements"""

print("=" * 80)
print("WATER BOX PERFORMANCE - ACTUAL MEASUREMENTS SUMMARY")
print("=" * 80)

# All measured data
data = [
    {'n': 8,    'ms': 0.003,  'steps_s': 335794, 'time_20k': 0.1},
    {'n': 27,   'ms': 0.031,  'steps_s': 31811,  'time_20k': 0.6},
    {'n': 64,   'ms': 0.159,  'steps_s': 6276,   'time_20k': 3.2},
    {'n': 125,  'ms': 0.643,  'steps_s': 1556,   'time_20k': 12.9},
    {'n': 216,  'ms': 1.751,  'steps_s': 571,    'time_20k': 35.0},
    {'n': 343,  'ms': 4.551,  'steps_s': 220,    'time_20k': 91.0},
    {'n': 512,  'ms': 9.730,  'steps_s': 103,    'time_20k': 194.6},
    {'n': 1000, 'ms': 28.248, 'steps_s': 35,     'time_20k': 565.0},
]

print("\n1. NON-DRUDE MODEL - ACTUAL PERFORMANCE")
print("-" * 50)
print(f"{'N waters':>8} | {'ms/step':>10} | {'steps/sec':>10} | {'20k steps':>12}")
print("-" * 50)

for d in data:
    if d['time_20k'] < 60:
        time_str = f"{d['time_20k']:.1f} sec"
    elif d['time_20k'] < 3600:
        time_str = f"{d['time_20k']/60:.1f} min"
    else:
        time_str = f"{d['time_20k']/3600:.1f} hr"
    
    print(f"{d['n']:8d} | {d['ms']:10.3f} | {d['steps_s']:10.0f} | {time_str:>12}")

print("\n\n2. DRUDE MODEL PREDICTIONS (Based on actual non-Drude times)")
print("-" * 70)
print("Assuming Drude is slower by different factors due to SCF iterations:")

for slowdown in [10, 20, 30, 50]:
    print(f"\n{slowdown}x slower:")
    print(f"{'N waters':>8} | {'ms/step':>10} | {'steps/sec':>10} | {'20k steps':>15}")
    print("-" * 60)
    
    for d in data:
        drude_ms = d['ms'] * slowdown
        drude_steps_s = 1000 / drude_ms
        drude_time_20k = drude_ms * 20000 / 1000
        
        if drude_time_20k < 60:
            time_str = f"{drude_time_20k:.1f} sec"
        elif drude_time_20k < 3600:
            time_str = f"{drude_time_20k/60:.1f} min"
        else:
            time_str = f"{drude_time_20k/3600:.1f} hr"
        
        print(f"{d['n']:8d} | {drude_ms:10.1f} | {drude_steps_s:10.1f} | {time_str:>15}")

print("\n\n3. PRACTICAL RECOMMENDATIONS")
print("-" * 50)
print("\nFor non-Drude simulations:")
print("- Up to 343 waters: < 2 minutes for 20k steps")
print("- Up to 512 waters: < 4 minutes for 20k steps")  
print("- 1000 waters: ~10 minutes for 20k steps")

print("\nFor Drude simulations (estimated):")
print("- 125 waters: 2-11 minutes (10-50x slower)")
print("- 216 waters: 6-29 minutes (10-50x slower)")
print("- 343 waters: 15-76 minutes (10-50x slower)")

print("\n\n4. ACTUAL SCALING BEHAVIOR")
print("-" * 50)
print("The scaling is slightly better than O(N²):")
print("- Small systems (8-64 waters): ~N^1.9")
print("- Medium systems (64-343 waters): ~N^1.95")
print("- Large systems (343-1000 waters): ~N^1.8")
print("\nThis is likely due to:")
print("- Cutoff radius limiting interactions")
print("- Cache effects for smaller systems")
print("- Memory bandwidth limitations for larger systems")

print("\n" + "=" * 80)