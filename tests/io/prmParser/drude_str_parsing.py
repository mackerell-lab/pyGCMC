#!/usr/bin/env python3
"""
测试PRM Parser是否能正确读取Drude STR文件中的力场参数
即使某些Drude特有参数无法解析也让它直接报错，后续再修正
"""

import pytest
import os
import sys

# Add pygcmc to path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../..'))

import pygcmc


@pytest.fixture
def str_file_path():
    """返回Drude STR文件路径"""
    return os.path.join(
        os.path.dirname(os.path.abspath(__file__)), 
        '../../../data/forcefields/charmm/c36_jul24/drude/drude_toppar_2023/toppar_drude_main_protein_2023a.str'
    )

def test_parse_drude_str_file(str_file_path):
    """测试是否能解析Drude STR文件"""
    # 尝试使用PRM Parser解析STR文件
    # 这里可能会报错，但没关系，让它报错
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # 如果能成功读取，检查基本信息
    assert ff is not None

def test_bond_parameters_from_str(str_file_path):
    """测试键参数的读取"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # 测试特定的键参数
    bond_params = ff.get_bond_params('ND2A2', 'CD2O1A')
    assert bond_params.kb == 376.20
    assert bond_params.b0 == 1.285
    
    bond_params = ff.get_bond_params('CD31C', 'CD32A')
    assert bond_params.kb == 222.50
    assert bond_params.b0 == 1.528  # 根据STR文件的实际值
    
    bond_params = ff.get_bond_params('ODW', 'HDW')
    assert bond_params.kb == 450.00
    assert bond_params.b0 == 0.9572

def test_angle_parameters_from_str(str_file_path):
    """测试角参数的读取"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # 测试特定的角参数
    # 使用实际在ANGLES部分第一个出现的参数
    angle_params = ff.get_angle_params('CD2O1A', 'ND2A2', 'CD31C')
    assert angle_params.ktheta == 40.90
    assert angle_params.theta0 == 116.10
    
    angle_params = ff.get_angle_params('HDP1A', 'ND2A1', 'HDP1A')
    assert angle_params.ktheta == 24.00
    assert angle_params.theta0 == 113.00

def test_dihedral_parameters_from_str(str_file_path):
    """测试二面角参数的读取"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # 测试特定的二面角参数
    dihedral_params = ff.get_dihedral_params('HDP1A', 'ND2A1', 'CD2O1A', 'OD2C1A')
    assert len(dihedral_params) > 0
    # 检查第一个参数
    assert dihedral_params[0].kchi == 2.000
    assert dihedral_params[0].n == 2
    assert dihedral_params[0].delta == 180.0

def test_improper_parameters_from_str(str_file_path):
    """测试不当二面角参数的读取"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # 测试特定的不当二面角参数
    # 使用实际存在的improper参数
    improper_params = ff.get_improper_params('CD2O1A', 'CD32A', 'ND2A2', 'OD2C1A')
    assert improper_params.kpsi == 100.00
    assert improper_params.psi0 == 0.0

def test_nonbonded_parameters_from_str(str_file_path):
    """测试非键参数的读取"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # 测试特定原子类型的LJ参数
    lj_params = ff.get_lj_params('HDP1A')
    assert lj_params.epsilon == -0.0100
    assert lj_params.rmin_half == 0.4000
    
    lj_params = ff.get_lj_params('CD31C')
    assert lj_params.epsilon == -0.0320
    assert lj_params.rmin_half == 1.8000
    
    lj_params = ff.get_lj_params('ND2A2')
    assert lj_params.epsilon == -0.2000
    assert lj_params.rmin_half == 1.8300
    
    lj_params = ff.get_lj_params('OD2C1A')
    assert lj_params.epsilon == -0.2000
    assert lj_params.rmin_half == 1.7800

def test_nbfix_parameters_from_str(str_file_path):
    """测试NBFIX参数的读取"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # 测试特定的NBFIX参数
    epsilon, rmin, found = ff.get_nbfix('ODW', 'CD2O3A')
    assert found == True
    assert epsilon == -0.11528  # 根据分析脚本的实际输出值
    assert rmin == 3.4869
    
    epsilon, rmin, found = ff.get_nbfix('ODW', 'ND2A2')
    assert found == True
    assert epsilon == -0.2054
    assert rmin == 3.6369

def test_drude_alpha_thole_parameters(str_file_path):
    """测试Drude ALPHA/THOLE参数的读取
    这个测试可能会失败，因为标准PRM parser可能不支持
    """
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # 尝试获取alpha和thole参数
    # 这些可能需要扩展ForceField类来支持
    # 让它直接报错，后续再修正
    alpha_params = ff.get_alpha_params('ODW')  # 可能不存在这个方法
    assert alpha_params.alpha == -0.97825258
    assert alpha_params.thole == 1.3
    
    alpha_params = ff.get_alpha_params('ND2A2')
    assert alpha_params.alpha == -1.858
    assert alpha_params.thole == 0.126

def test_lonepair_definitions(str_file_path):
    """测试LONEPAIR定义的读取
    这个测试可能会失败，因为标准PRM parser可能不支持
    """
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # 尝试获取lonepair定义
    # 这些可能需要扩展ForceField类来支持
    lonepairs = ff.get_lonepairs()  # 可能不存在这个方法
    assert len(lonepairs) == 143

def test_anisotropy_definitions(str_file_path):
    """测试ANISOTROPY定义的读取
    这个测试可能会失败，因为标准PRM parser可能不支持
    """
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # 尝试获取anisotropy定义
    # 这些可能需要扩展ForceField类来支持
    anisotropies = ff.get_anisotropies()  # 可能不存在这个方法
    assert len(anisotropies) == 68

def test_atom_type_count(str_file_path):
    """测试原子类型数量"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # 检查是否读取了所有原子类型
    # 根据分析结果，应该有175个原子类型
    # 使用get_num_lj_params来计算原子类型数量
    num_atom_types = ff.get_num_lj_params()
    assert num_atom_types >= 175

def test_specific_drude_atom_types(str_file_path):
    """测试特定的Drude原子类型"""
    ff = pygcmc.PRMParser.parse_file(str_file_path)
    
    # 测试Drude粒子类型
    lj_params = ff.get_lj_params('DRUD')
    assert lj_params.epsilon == 0.0
    assert lj_params.rmin_half == 0.01
    
    # 测试水的Drude粒子
    lj_params = ff.get_lj_params('DOH2')
    assert lj_params.epsilon == 0.0
    assert lj_params.rmin_half == 0.01
    
    # 测试孤对电子
    lj_params = ff.get_lj_params('LPD')
    assert lj_params.epsilon == 0.0
    assert lj_params.rmin_half == 0.01
    
    lj_params = ff.get_lj_params('LPDW')
    assert lj_params.epsilon == 0.0
    assert lj_params.rmin_half == 0.01

if __name__ == "__main__":
    pytest.main([__file__, "-v"])