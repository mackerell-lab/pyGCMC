"""
测试装饰器，用于处理全局状态重置和其他测试配置
"""

import functools
import pytest
import pygcmc
import gc


def reset_pgp_state(func):
    """装饰器：在测试前后重置PGP状态"""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # 测试前重置
        if hasattr(pygcmc, 'resetPGPState'):
            pygcmc.resetPGPState()
        if hasattr(pygcmc, 'clearPMEState'):
            pygcmc.clearPMEState()
        
        try:
            # 运行测试
            result = func(*args, **kwargs)
        finally:
            # 测试后重置
            if hasattr(pygcmc, 'resetPGPState'):
                pygcmc.resetPGPState()
            if hasattr(pygcmc, 'clearPMEState'):
                pygcmc.clearPMEState()
            gc.collect()
        
        return result
    return wrapper


def serial_execution(func):
    """标记测试为串行执行"""
    return pytest.mark.serial(func)


def pgp_test(func):
    """组合装饰器：用于PGP测试"""
    @functools.wraps(func)
    @reset_pgp_state
    @pytest.mark.xdist_group("pgp")
    def wrapper(*args, **kwargs):
        return func(*args, **kwargs)
    return wrapper


def unstable_test(reason="Known memory issue"):
    """标记不稳定的测试"""
    def decorator(func):
        @functools.wraps(func)
        @reset_pgp_state
        @pytest.mark.xdist_group("unstable")
        @pytest.mark.flaky(reruns=2, reruns_delay=1)
        def wrapper(*args, **kwargs):
            print(f"\nRunning unstable test: {reason}")
            return func(*args, **kwargs)
        return wrapper
    return decorator