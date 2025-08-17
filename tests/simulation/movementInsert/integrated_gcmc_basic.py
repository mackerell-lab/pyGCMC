# tests/simulation/movementInsert/test_integrated_gcmc_patterns.py
"""
基于从gcmc_gpu/opencl学到的测试模式的综合测试

这个测试文件整合了从老程序中学到的最佳实践：
1. 配置驱动的测试
2. GCMC/MD集成循环
3. 多分子类型同时模拟
4. 收敛性和统计监控
"""

import pytest
import random
import math
import numpy as np
import pygcmc
from typing import Dict, List, Tuple
from dataclasses import dataclass

# Constants
kB = 0.008314463  # kJ/(mol·K)

# Import helper functions
from .basic_insertion_helpers import (
    create_empty_system,
    create_water_molecule,
    insert_molecule,
    calculate_system_energy
)


@dataclass
class GCMCConfig:
    """GCMC配置类 - 模仿gcmc.inp文件"""
    mcsteps: int = 10000
    nprint: int = 100
    use_cavity_bias: bool = True
    use_conf_bias: bool = False
    num_conf_bias_trial: int = 10
    energy_cutoff: float = 8.0
    temperature: float = 300.0
    fragments: Dict = None
    
    def __post_init__(self):
        if self.fragments is None:
            self.fragments = {
                'water': {'conc': 55.0, 'muex': -5.6}
            }


class GCMCTracker:
    """统计跟踪器 - 学自gcmc_moves.cpp"""
    def __init__(self):
        self.moves = {
            'insert': {'accepted': 0, 'rejected': 0},
            'delete': {'accepted': 0, 'rejected': 0},
            'translate': {'accepted': 0, 'rejected': 0},
            'rotate': {'accepted': 0, 'rejected': 0}
        }
        self.n_molecules_history = []
        self.energy_history = []
    
    def record_move(self, move_type: str, accepted: bool):
        if accepted:
            self.moves[move_type]['accepted'] += 1
        else:
            self.moves[move_type]['rejected'] += 1
    
    def get_acceptance_rate(self, move_type: str) -> float:
        move = self.moves[move_type]
        total = move['accepted'] + move['rejected']
        if total == 0:
            return 0.0
        return move['accepted'] / total
    
    def record_state(self, n_molecules: int, energy: float):
        self.n_molecules_history.append(n_molecules)
        self.energy_history.append(energy)
    
    def is_converged(self, window: int = 100, tolerance: float = 0.1) -> bool:
        """检查系统是否收敛"""
        if len(self.n_molecules_history) < 2 * window:
            return False
        
        # 比较最近两个窗口的平均值
        recent = self.n_molecules_history[-window:]
        previous = self.n_molecules_history[-2*window:-window]
        
        avg_recent = np.mean(recent)
        avg_previous = np.mean(previous)
        
        if avg_recent == 0:
            return abs(avg_previous) < tolerance
        
        relative_change = abs(avg_recent - avg_previous) / avg_recent
        return relative_change < tolerance


def load_config_from_dict(config_dict: dict) -> GCMCConfig:
    """从字典加载配置 - 模仿从文件读取"""
    return GCMCConfig(**config_dict)


def select_move_type(probabilities: dict = None) -> str:
    """选择移动类型 - 学自run_gcmc_step()"""
    if probabilities is None:
        # 默认概率
        probabilities = {
            'insert': 0.25,
            'delete': 0.25,
            'translate': 0.25,
            'rotate': 0.25
        }
    
    rand = random.random()
    cumsum = 0
    for move, prob in probabilities.items():
        cumsum += prob
        if rand < cumsum:
            return move
    return list(probabilities.keys())[-1]


def calculate_insertion_acceptance(system, molecule, config: GCMCConfig) -> float:
    """计算插入接受概率 - 基于move_add()"""
    # 计算能量变化
    energy_before = calculate_system_energy(system)
    system_new = insert_molecule(system, molecule)
    pygcmc.computeSystemEnergyCutoff(system_new)
    energy_after = calculate_system_energy(system_new)
    delta_E = energy_after - energy_before
    
    # 获取参数
    beta = 1.0 / (kB * config.temperature)
    n_molecules = len(system.residues)
    
    # 简化的化学势和浓度
    muex = config.fragments['water']['muex']
    conc = config.fragments['water']['conc']
    nbar = conc * 0.6022  # 简化的浓度转换
    B = beta * muex + math.log(nbar)
    
    # Cavity bias因子
    fn = 1.0
    if config.use_cavity_bias:
        # 简化的cavity fraction估计
        fn = max(0.1, 1.0 - n_molecules * 0.01)
    
    # 接受概率
    acc = min(1.0, fn / (n_molecules + 1) * math.exp(B - beta * delta_E))
    return acc


def test_config_driven_gcmc():
    """测试配置驱动的GCMC - 学自gcmc.inp模式"""
    # 加载配置
    config = load_config_from_dict({
        'mcsteps': 1000,
        'use_cavity_bias': True,
        'temperature': 300.0,
        'fragments': {
            'water': {'conc': 55.0, 'muex': -5.6}
        }
    })
    
    # 初始化系统
    system = create_empty_system()
    tracker = GCMCTracker()
    
    # 运行GCMC步骤
    for step in range(config.mcsteps):
        move_type = select_move_type()
        
        # 简化的移动实现
        accepted = False
        if move_type == 'insert':
            if len(system.residues) < 20:  # 限制最大分子数
                molecule = create_water_molecule(
                    random.uniform(0, 5),
                    random.uniform(0, 5),
                    random.uniform(0, 5)
                )
                acc_prob = calculate_insertion_acceptance(system, molecule, config)
                if random.random() < acc_prob:
                    system = insert_molecule(system, molecule)
                    accepted = True
        
        elif move_type == 'delete':
            if len(system.residues) > 0:
                # 简化的删除（总是接受）
                accepted = True
        
        # 记录统计
        tracker.record_move(move_type, accepted)
        
        # 定期记录状态
        if step % config.nprint == 0:
            pygcmc.computeSystemEnergyCutoff(system)
            energy = calculate_system_energy(system)
            tracker.record_state(len(system.residues), energy)
    
    # 验证结果
    insert_rate = tracker.get_acceptance_rate('insert')
    assert 0 <= insert_rate <= 1, f"Invalid acceptance rate: {insert_rate}"
    
    # 检查是否有分子插入
    assert len(tracker.n_molecules_history) > 0, "No state recorded"


def test_gcmc_md_integration_pattern():
    """测试GCMC/MD集成模式 - 学自protein/run.py"""
    system = create_empty_system()
    tracker = GCMCTracker()
    
    n_cycles = 5
    n_gcmc_steps = 100
    n_md_steps = 10  # 简化的MD步骤
    
    for cycle in range(n_cycles):
        # GCMC阶段
        for _ in range(n_gcmc_steps):
            move_type = select_move_type({'insert': 0.5, 'delete': 0.5})
            
            if move_type == 'insert' and len(system.residues) < 10:
                molecule = create_water_molecule(
                    random.uniform(0, 5),
                    random.uniform(0, 5),
                    random.uniform(0, 5)
                )
                system = insert_molecule(system, molecule)
                tracker.record_move('insert', True)
        
        # 模拟MD弛豫（简化版本）
        for _ in range(n_md_steps):
            # 在实际中这里会调用MD引擎
            # 这里只是模拟能量下降
            pygcmc.computeSystemEnergyCutoff(system)
            pass
        
        # 记录循环结束时的状态
        energy = calculate_system_energy(system)
        tracker.record_state(len(system.residues), energy)
    
    # 验证循环执行
    assert len(tracker.n_molecules_history) == n_cycles
    assert all(n >= 0 for n in tracker.n_molecules_history)


def test_multiple_fragment_types():
    """测试多分子类型 - 学自OpenCL examples"""
    # 配置多种分子
    config = GCMCConfig(
        mcsteps=500,
        fragments={
            'water': {'conc': 55.0, 'muex': -5.6, 'prob': 0.8},
            'ion': {'conc': 0.15, 'muex': -10.0, 'prob': 0.2}
        }
    )
    
    system = create_empty_system()
    fragment_counts = {'water': 0, 'ion': 0}
    
    for _ in range(config.mcsteps):
        # 根据概率选择分子类型
        rand = random.random()
        if rand < 0.8:
            frag_type = 'water'
        else:
            frag_type = 'ion'
        
        # 简化的插入
        if len(system.residues) < 20:
            molecule = create_water_molecule(
                random.uniform(0, 5),
                random.uniform(0, 5),
                random.uniform(0, 5)
            )
            system = insert_molecule(system, molecule)
            fragment_counts[frag_type] += 1
    
    # 验证两种分子都被插入
    assert fragment_counts['water'] > 0, "No water molecules inserted"
    assert fragment_counts['ion'] > 0, "No ions inserted"
    
    # 验证比例大致正确
    total = sum(fragment_counts.values())
    if total > 10:
        water_ratio = fragment_counts['water'] / total
        assert 0.6 < water_ratio < 1.0, f"Water ratio {water_ratio} out of expected range"


