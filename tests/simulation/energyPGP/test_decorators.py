"""
PGP测试专用装饰器
"""

import functools
import pytest
import pygcmc
import gc
import os


def force_serial(func):
    """
    装饰器：强制测试串行执行
    使用 pytest.mark.serial 标记，配合 pytest-xdist 的分组功能
    """
    @functools.wraps(func)
    @pytest.mark.serial  # 标记为串行执行
    @pytest.mark.xdist_group("serial_only")  # 确保在同一个worker中执行
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper


def reset_state(func):
    """
    装饰器：在测试前后彻底重置状态
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # 测试前重置
        try:
            if hasattr(pygcmc, 'resetPGPState'):
                pygcmc.resetPGPState()
            if hasattr(pygcmc, 'clearPMEState'):
                pygcmc.clearPMEState()
        except:
            pass
        
        gc.collect()
        
        try:
            result = func(*args, **kwargs)
        finally:
            # 测试后重置
            try:
                if hasattr(pygcmc, 'resetPGPState'):
                    pygcmc.resetPGPState()
                if hasattr(pygcmc, 'clearPMEState'):
                    pygcmc.clearPMEState()
            except:
                pass
            gc.collect()
        
        return result
    return wrapper


def pgp_unstable_test(func):
    """
    组合装饰器：用于不稳定的PGP测试
    通过多次重置和垃圾回收来减少内存问题
    """
    @functools.wraps(func)
    @reset_state
    def wrapper(*args, **kwargs):
        print(f"\n[特殊处理] {func.__name__}")
        
        # 在执行前进行额外的清理
        for _ in range(2):
            if hasattr(pygcmc, 'resetPGPState'):
                pygcmc.resetPGPState()
            gc.collect()
        
        try:
            result = func(*args, **kwargs)
        finally:
            # 执行后立即清理
            if hasattr(pygcmc, 'resetPGPState'):
                pygcmc.resetPGPState()
            gc.collect()
            
        return result
    return wrapper