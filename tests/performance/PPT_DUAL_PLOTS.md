# GCMC Performance Comparison: Additive vs Polarizable Water Models

## Performance Benchmark (steps/second)

| Waters | TIP3P (Additive) | SWM4-NDP (Polarizable) | Speed Ratio |
|:------:|:----------------:|:----------------------:|:-----------:|
| 16     | 55,246           | 3,267                  | 16.9x       |
| 32     | 11,976           | 529                    | 22.6x       |
| 64     | 3,833            | 90                     | 42.7x       |
| 128    | 734              | 14                     | 53.0x       |
| 256    | 205              | 2                      | 99.7x       |

## Scaling Analysis - Two Views

### Figure 1: Log-Log Scaling
[Insert Log-Log Plot Here]

### Figure 2: N² Scaling Analysis  
[Insert N² Plot Here]

```python
import matplotlib.pyplot as plt
import numpy as np

# Data
N = np.array([16, 32, 64, 128, 256])
TIP3P = np.array([55246, 11976, 3833, 734, 205])
SWM4 = np.array([3267, 529, 90, 14, 2])

# Convert to time per step (milliseconds)
time_TIP3P = 1000.0 / TIP3P  # ms per step
time_SWM4 = 1000.0 / SWM4    # ms per step

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# ====== Plot 1: Log-Log ======
ax1.loglog(N, TIP3P, 'o-', label='TIP3P', markersize=10, linewidth=2.5, color='#2E86AB')
ax1.loglog(N, SWM4, 's-', label='SWM4-NDP', markersize=10, linewidth=2.5, color='#A23B72')

# Power law fits
p_tip3p = np.polyfit(np.log(N), np.log(TIP3P), 1)
p_swm4 = np.polyfit(np.log(N), np.log(SWM4), 1)

# Plot fit lines
N_fit = np.logspace(1.2, 2.4, 100)
ax1.loglog(N_fit, np.exp(p_tip3p[1]) * N_fit**p_tip3p[0], 
           '--', color='#2E86AB', alpha=0.5, 
           label=f'TIP3P: $N^{{{p_tip3p[0]:.2f}}}$')
ax1.loglog(N_fit, np.exp(p_swm4[1]) * N_fit**p_swm4[0], 
           '--', color='#A23B72', alpha=0.5,
           label=f'SWM4: $N^{{{p_swm4[0]:.2f}}}$')

ax1.set_xlabel('Number of Waters (N)', fontsize=14)
ax1.set_ylabel('GCMC Steps/Second', fontsize=14)
ax1.set_title('(a) Log-Log Scaling', fontsize=16)
ax1.legend(fontsize=12, loc='upper right')
ax1.grid(True, alpha=0.3)

# Add speed ratio annotations
for i, (n, ratio) in enumerate(zip(N, TIP3P/SWM4)):
    ax1.annotate(f'{ratio:.0f}x', 
                xy=(n, np.sqrt(TIP3P[i]*SWM4[i])), 
                fontsize=10, ha='center', va='bottom')

# ====== Plot 2: N² Scaling ======
N_squared = N**2

# Plot time per step vs N²
ax2.plot(N_squared, time_TIP3P, 'o-', label='TIP3P', 
         markersize=10, linewidth=2.5, color='#2E86AB')
ax2.plot(N_squared, time_SWM4, 's-', label='SWM4-NDP', 
         markersize=10, linewidth=2.5, color='#A23B72')

# Fit lines for ideal O(N²) scaling
# For perfect O(N²), time = a*N²
fit_tip3p = np.polyfit(N_squared, time_TIP3P, 1)
fit_swm4_linear = np.polyfit(N_squared[:3], time_SWM4[:3], 1)  # Fit only small systems

# Plot ideal O(N²) lines
N_squared_fit = np.linspace(0, 70000, 100)
ax2.plot(N_squared_fit, fit_tip3p[0] * N_squared_fit + fit_tip3p[1], 
         '--', color='#2E86AB', alpha=0.5, label='TIP3P linear fit')
ax2.plot(N_squared_fit, fit_swm4_linear[0] * N_squared_fit + fit_swm4_linear[1], 
         ':', color='#A23B72', alpha=0.5, label='SWM4 O(N²) reference')

ax2.set_xlabel('N² (Number of Waters Squared)', fontsize=14)
ax2.set_ylabel('Time per GCMC Step (ms)', fontsize=14)
ax2.set_title('(b) Deviation from O(N²) Scaling', fontsize=16)
ax2.legend(fontsize=12)
ax2.grid(True, alpha=0.3)

# Add annotations showing deviation
ax2.text(40000, 100, 'SWM4 deviates from\nO(N²) for large N', 
         fontsize=11, ha='center', bbox=dict(boxstyle="round,pad=0.3", 
         facecolor='yellow', alpha=0.5))

# Set y-axis to start from 0
ax2.set_ylim(bottom=0)
ax2.set_xlim(left=0)

plt.tight_layout()
plt.savefig('gcmc_performance_dual_analysis.png', dpi=300, bbox_inches='tight')
plt.show()

# ====== Additional Analysis Plot: Time Ratio vs N ======
plt.figure(figsize=(8, 6))
time_ratio = time_SWM4 / time_TIP3P
plt.plot(N, time_ratio, 'o-', markersize=10, linewidth=2.5, color='#8B4513')

# Fit exponential growth
p_ratio = np.polyfit(np.log(N), np.log(time_ratio), 1)
N_fit = np.linspace(16, 256, 100)
plt.plot(N_fit, np.exp(p_ratio[1]) * N_fit**p_ratio[0], 
         '--', color='#8B4513', alpha=0.5,
         label=f'Ratio ∝ $N^{{{p_ratio[0]:.2f}}}$')

plt.xlabel('Number of Waters (N)', fontsize=14)
plt.ylabel('Time Ratio (SWM4/TIP3P)', fontsize=14)
plt.title('Performance Gap Growth with System Size', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)

# Add values
for n, ratio in zip(N, time_ratio):
    plt.annotate(f'{ratio:.0f}x', xy=(n, ratio), 
                xytext=(5, 5), textcoords='offset points')

plt.tight_layout()
plt.savefig('gcmc_performance_ratio_growth.png', dpi=300, bbox_inches='tight')
plt.show()
```

## Key Scientific Insights

### From Log-Log Plot:
• Both models show power-law scaling
• TIP3P: Nearly perfect O(N²·⁰²)
• SWM4: Significant deviation at O(N²·⁶⁵)

### From N² Plot:
• **TIP3P**: Nearly linear relationship with N² (ideal O(N²) scaling)
• **SWM4**: Strong upward curvature, indicating super-quadratic scaling
• Deviation becomes dramatic for N > 100

### Performance Gap Analysis:
• Speed ratio grows as ~N⁰·⁶³
• Extrapolation: 1000 waters → ~200x slower
• SWM4 becomes impractical for production GCMC

---

*Test Environment: Single-core CPU | PyGCMC | 1.2 nm cutoff | C++ implementation*