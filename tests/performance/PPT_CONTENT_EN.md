# GCMC Performance Comparison: Additive vs Polarizable Water Models

## Performance Benchmark (steps/second)

| Waters | TIP3P (Additive) | SWM4-NDP (Polarizable) | Speed Ratio |
|:------:|:----------------:|:----------------------:|:-----------:|
| 8      | 389,805          | 7,426                  | 52.5x       |
| 27     | 21,136           | 799                    | 26.5x       |
| 64     | 3,722            | 84                     | 44.3x       |
| 125    | 1,041            | 13                     | 77.4x       |

**Average Speed Difference: 50x**

## Scaling Behavior

[Insert Figure Here]

```python
# Plotting code
import matplotlib.pyplot as plt
import numpy as np

# Data
N = np.array([8, 27, 64, 125])
TIP3P = np.array([389805, 21136, 3722, 1041])
SWM4 = np.array([7426, 799, 84, 13])

# Create log-log plot
plt.figure(figsize=(8, 6))
plt.loglog(N, TIP3P, 'o-', label='TIP3P', markersize=10, linewidth=2, color='#2E86AB')
plt.loglog(N, SWM4, 's-', label='SWM4-NDP', markersize=10, linewidth=2, color='#A23B72')

# Add power law fits
p_tip3p = np.polyfit(np.log(N), np.log(TIP3P), 1)
p_swm4 = np.polyfit(np.log(N), np.log(SWM4), 1)

# Plot fit lines
N_fit = np.logspace(0.9, 2.1, 100)
plt.loglog(N_fit, np.exp(p_tip3p[1]) * N_fit**p_tip3p[0], 
           '--', color='#2E86AB', alpha=0.5, 
           label=f'TIP3P: N^{-p_tip3p[0]:.1f}')
plt.loglog(N_fit, np.exp(p_swm4[1]) * N_fit**p_swm4[0], 
           '--', color='#A23B72', alpha=0.5,
           label=f'SWM4: N^{-p_swm4[0]:.1f}')

plt.xlabel('Number of Waters (N)', fontsize=14)
plt.ylabel('GCMC Steps/Second', fontsize=14)
plt.title('GCMC Performance Scaling', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)
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
| Scaling | O(N²·²) | O(N²·³) |

## Key Findings

• **SWM4 is 26-77x slower than TIP3P**
• **Both scale approximately as O(N²)**
• **64 water system: 3,722 vs 84 steps/sec**
• **1M steps: 5 minutes vs 3.3 hours**

## Recommendations

- **Small systems (<100 waters)**: Both viable, SWM4 for higher accuracy
- **Medium systems (100-500 waters)**: Prefer TIP3P
- **Large systems (>500 waters)**: TIP3P only or consider GPU acceleration

---

*Test Environment: Single-core CPU | PyGCMC | 1.2 nm cutoff | C++ implementation*