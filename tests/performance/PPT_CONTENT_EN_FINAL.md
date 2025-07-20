# GCMC Performance Comparison: Additive vs Polarizable Water Models

## Performance Benchmark (steps/second)

| Waters | TIP3P (Additive) | SWM4-NDP (Polarizable) | Speed Ratio |
|:------:|:----------------:|:----------------------:|:-----------:|
| 16     | 55,246           | 3,267                  | 16.9x       |
| 32     | 11,976           | 529                    | 22.6x       |
| 64     | 3,833            | 90                     | 42.7x       |
| 128    | 734              | 14                     | 53.0x       |
| 256    | 205              | 2                      | 99.7x       |

**Speed ratio increases monotonically: approaching 100x for large systems**

## Scaling Behavior

[Insert Figure Here]

```python
# Plotting code
import matplotlib.pyplot as plt
import numpy as np

# Data
N = np.array([16, 32, 64, 128, 256])
TIP3P = np.array([55246, 11976, 3833, 734, 205])
SWM4 = np.array([3267, 529, 90, 14, 2])

# Create log-log plot
plt.figure(figsize=(9, 7))
plt.loglog(N, TIP3P, 'o-', label='TIP3P', markersize=10, linewidth=2.5, color='#2E86AB')
plt.loglog(N, SWM4, 's-', label='SWM4-NDP', markersize=10, linewidth=2.5, color='#A23B72')

# Add power law fits
p_tip3p = np.polyfit(np.log(N), np.log(TIP3P), 1)
p_swm4 = np.polyfit(np.log(N), np.log(SWM4), 1)

# Plot fit lines
N_fit = np.logspace(1.2, 2.4, 100)
plt.loglog(N_fit, np.exp(p_tip3p[1]) * N_fit**p_tip3p[0], 
           '--', color='#2E86AB', alpha=0.5, 
           label=f'TIP3P: $N^{{{p_tip3p[0]:.2f}}}$')
plt.loglog(N_fit, np.exp(p_swm4[1]) * N_fit**p_swm4[0], 
           '--', color='#A23B72', alpha=0.5,
           label=f'SWM4: $N^{{{p_swm4[0]:.2f}}}$')

plt.xlabel('Number of Waters (N)', fontsize=14)
plt.ylabel('GCMC Steps/Second', fontsize=14)
plt.title('GCMC Performance Scaling', fontsize=16)
plt.legend(fontsize=12, loc='upper right')
plt.grid(True, alpha=0.3)

# Add annotations for speed ratios
for i, (n, ratio) in enumerate(zip(N, TIP3P/SWM4)):
    plt.annotate(f'{ratio:.0f}x', 
                xy=(n, np.sqrt(TIP3P[i]*SWM4[i])), 
                fontsize=10, ha='center', va='bottom')

# Highlight 100x threshold
plt.axhline(y=2, color='red', linestyle=':', alpha=0.5)
plt.text(20, 2.5, 'SWM4 @ 2 steps/s', fontsize=9, color='red')

plt.tight_layout()
plt.savefig('gcmc_performance_scaling.png', dpi=300, bbox_inches='tight')
plt.show()
```

## Model Comparison

| Feature | TIP3P | SWM4-NDP |
|:-------:|:-----:|:--------:|
| Sites/water | 3 | 5 |
| Polarization | None | SCF iteration |
| Accuracy | ++ | ++++ |
| Scaling | O(N²·⁰²) | O(N²·⁶⁵) |

## Key Findings

• **SWM4 is 17-100x slower than TIP3P**
• **Speed gap widens dramatically with system size**
• **256 waters: ~100x performance difference**
• **SWM4 scales poorly: O(N²·⁶⁵) vs O(N²·⁰²)**
• **Large systems: SWM4 becomes impractical (2 steps/sec)**

## Recommendations

- **Small systems (<50 waters)**: Both viable, ~20x difference
- **Medium systems (50-200 waters)**: TIP3P preferred, 40-50x faster
- **Large systems (>200 waters)**: TIP3P only, approaching 100x faster
- **GCMC applications**: TIP3P strongly recommended due to frequent energy calculations

---

*Test Environment: Single-core CPU | PyGCMC | 1.2 nm cutoff | C++ implementation*  
*Note: SWM4 SCF convergence becomes challenging for N>200*