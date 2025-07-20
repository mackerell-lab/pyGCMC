# GCMC性能对比：加性vs极化水模型

## 性能对比 (步/秒)

| 水分子数 | TIP3P (加性) | SWM4-NDP (极化) | 速度比 |
|:--------:|:------------:|:---------------:|:------:|
| 8        | 389,805      | 7,426           | 52.5x  |
| 27       | 21,136       | 799             | 26.5x  |
| 64       | 3,722        | 84              | 44.3x  |
| 125      | 1,041        | 13              | 77.4x  |

**平均速度差异: 50倍**

## 缩放曲线

[此处插入图表]

```python
# 作图代码
import matplotlib.pyplot as plt
import numpy as np

# 数据
N = np.array([8, 27, 64, 125])
TIP3P = np.array([389805, 21136, 3722, 1041])
SWM4 = np.array([7426, 799, 84, 13])

# 创建log-log图
plt.figure(figsize=(8, 6))
plt.loglog(N, TIP3P, 'o-', label='TIP3P', markersize=10, linewidth=2, color='#2E86AB')
plt.loglog(N, SWM4, 's-', label='SWM4-NDP', markersize=10, linewidth=2, color='#A23B72')

# 添加幂律拟合
p_tip3p = np.polyfit(np.log(N), np.log(TIP3P), 1)
p_swm4 = np.polyfit(np.log(N), np.log(SWM4), 1)

# 绘制拟合线
N_fit = np.logspace(0.9, 2.1, 100)
plt.loglog(N_fit, np.exp(p_tip3p[1]) * N_fit**p_tip3p[0], 
           '--', color='#2E86AB', alpha=0.5, 
           label=f'TIP3P: N^{-p_tip3p[0]:.1f}')
plt.loglog(N_fit, np.exp(p_swm4[1]) * N_fit**p_swm4[0], 
           '--', color='#A23B72', alpha=0.5,
           label=f'SWM4: N^{-p_swm4[0]:.1f}')

plt.xlabel('水分子数 (N)', fontsize=14)
plt.ylabel('GCMC步数/秒', fontsize=14)
plt.title('GCMC性能缩放', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('gcmc_performance_scaling.png', dpi=300, bbox_inches='tight')
plt.show()
```

## 模型特性对比

| 特性 | TIP3P | SWM4-NDP |
|:----:|:-----:|:--------:|
| 原子/水 | 3 | 5 |
| 极化 | 无 | SCF迭代 |
| 精度 | ++ | ++++ |
| 缩放 | O(N²·²) | O(N²·³) |

## 关键发现

• **SWM4比TIP3P慢26-77倍**
• **两者都接近O(N²)缩放**
• **64水系统：3722 vs 84步/秒**
• **100万步：5分钟 vs 3.3小时**

## 应用建议

- **小系统(<100水)**: 两者皆可，SWM4提供更高精度
- **中等系统(100-500水)**: 优先选择TIP3P
- **大系统(>500水)**: 仅使用TIP3P或考虑GPU加速

---

*测试环境: 单核CPU | PyGCMC | 截止1.2nm | C++实现*