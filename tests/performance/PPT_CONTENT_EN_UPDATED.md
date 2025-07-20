# GCMC Performance Comparison: Additive vs Polarizable Water Models

## Performance Benchmark (steps/second)

| Waters | TIP3P (Additive) | SWM4-NDP (Polarizable) | Speed Ratio |
|:------:|:----------------:|:----------------------:|:-----------:|
| 16     | 55,246           | 3,267                  | 16.9x       |
| 32     | 11,976           | 529                    | 22.6x       |
| 64     | 3,833            | 90                     | 42.7x       |
| 128    | 734              | 14                     | 53.0x       |

**Speed ratio increases monotonically with system size**

## Scaling Behavior

[Insert Figure Here]

```python
# Plotting code
import matplotlib.pyplot as plt
import numpy as np

# Data
N = np.array([16, 32, 64, 128])
TIP3P = np.array([55246, 11976, 3833, 734])
SWM4 = np.array([3267, 529, 90, 14])

# Create log-log plot
plt.figure(figsize=(8, 6))
plt.loglog(N, TIP3P, 'o-', label='TIP3P', markersize=10, linewidth=2, color='#2E86AB')
plt.loglog(N, SWM4, 's-', label='SWM4-NDP', markersize=10, linewidth=2, color='#A23B72')

# Add power law fits
p_tip3p = np.polyfit(np.log(N), np.log(TIP3P), 1)
p_swm4 = np.polyfit(np.log(N), np.log(SWM4), 1)

# Plot fit lines
N_fit = np.logspace(1.2, 2.1, 100)
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

• **SWM4 is 17-53x slower than TIP3P**
• **Speed gap widens with system size**
• **TIP3P scales as ~O(N²)**
• **SWM4 scales as ~O(N²·⁶⁵)** due to SCF overhead
• **128 waters: 734 vs 14 steps/sec**

## Recommendations

- **Small systems (<50 waters)**: Both viable, SWM4 for accuracy
- **Medium systems (50-200 waters)**: Strongly prefer TIP3P  
- **Large systems (>200 waters)**: TIP3P only (100x faster)
- **GCMC applications**: TIP3P strongly recommended

---

*Test Environment: Single-core CPU | PyGCMC | 1.2 nm cutoff | C++ implementation*