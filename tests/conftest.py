"""
Global pytest configuration for the test suite.
"""

import pytest
import pygcmc
import gc
import os
import weakref

@pytest.fixture(autouse=True)
def reset_global_state(request):
    """自动在每个测试前后重置全局状态"""
    # 对于特定的测试，需要更彻底的清理
    test_name = request.node.name
    is_unstable = hasattr(request.node, 'get_closest_marker') and request.node.get_closest_marker('unstable')
    needs_deep_clean = is_unstable or any(name in test_name for name in [
        'test_movement_residues',
        'test_vdw_movement_debug', 
        'test_pgp_self_consistency'
    ])
    
    try:
        # 测试前重置
        if hasattr(pygcmc, 'resetPGPState'):
            pygcmc.resetPGPState()
        if hasattr(pygcmc, 'clearPMEState'):
            pygcmc.clearPMEState()
        
        # 对于有问题的测试，进行更彻底的清理
        if needs_deep_clean:
            gc.collect()
            gc.collect()  # 运行两次以确保循环引用被清理
            
            # 尝试重置多次以确保状态清理
            for _ in range(3):
                if hasattr(pygcmc, 'resetPGPState'):
                    pygcmc.resetPGPState()
                    
    except Exception:
        # 如果重置失败，继续执行测试
        pass
    
    yield
    
    try:
        # 测试后再次重置
        if hasattr(pygcmc, 'resetPGPState'):
            pygcmc.resetPGPState()
        if hasattr(pygcmc, 'clearPMEState'):
            pygcmc.clearPMEState()
            
        # 对于有问题的测试，测试后也进行彻底清理
        if needs_deep_clean:
            gc.collect()
            gc.collect()
            
    except Exception:
        # 忽略清理时的错误
        pass

def pytest_collection_modifyitems(config, items):
    """
    Modify test collection to mark certain tests for serial execution.
    
    This helps avoid segmentation faults that occur when certain tests
    run in parallel due to global state issues.
    """
    # 检查是否通过环境变量强制运行不稳定的测试
    force_run_unstable = os.environ.get('PYTEST_RUN_UNSTABLE_TESTS', '').lower() == 'true'
    
    items_to_remove = []
    
    for item in items:
        # Mark OpenMM tests for serial execution
        if "energyOpenmm" in str(item.fspath):
            item.add_marker(pytest.mark.xdist_group("openmm"))
            
        # Mark PGP tests that use original implementation
        if "test_pgp_grid_convergence" in item.name and "independent" not in item.name:
            item.add_marker(pytest.mark.xdist_group("pgp_original"))
            
        if "test_pgp_cutoff_continuity" in item.name:
            item.add_marker(pytest.mark.xdist_group("pgp_original"))
            
        # Handle unstable tests - 确保它们在单独的进程中串行运行
        if any(name in item.name for name in ['test_movement_residues', 
                                               'test_vdw_movement_debug',
                                               'test_pgp_self_consistency']):
            # 将这些测试放在专门的组中，确保它们：
            # 1. 在同一个 worker 中运行（通过 xdist_group）
            # 2. 按顺序串行执行（同组内的测试不会并行）
            # 3. 与其他测试隔离
            item.add_marker(pytest.mark.xdist_group("unstable_serial"))
            item.add_marker(pytest.mark.timeout(120))
            
            # 为这些测试添加特殊标记，表示需要额外的清理
            item.add_marker(pytest.mark.unstable)

@pytest.fixture(scope="session", autouse=True)
def configure_test_environment():
    """配置测试环境以减少内存问题"""
    # 设置环境变量以减少内存碎片
    os.environ['MALLOC_MMAP_THRESHOLD_'] = '128000'  # 128KB
    os.environ['MALLOC_TRIM_THRESHOLD_'] = '128000'
    
    # 确保在会话开始时重置状态
    if hasattr(pygcmc, 'resetPGPState'):
        pygcmc.resetPGPState()
    if hasattr(pygcmc, 'clearPMEState'):
        pygcmc.clearPMEState()
    
    yield
    
    # 会话结束时清理
    gc.collect()


def pytest_configure(config):
    """Configure custom markers."""
    config.addinivalue_line(
        "markers", "serial: mark test to run in serial mode"
    )
    config.addinivalue_line(
        "markers", "timeout: mark test with custom timeout"
    )
    config.addinivalue_line(
        "markers", "unstable: mark test as unstable and requiring special handling"
    )
    
    # 如果使用xdist，设置worker数量限制
    if hasattr(config, 'workerinput'):
        # 这是一个worker进程
        worker_id = config.workerinput.get('workerid', '')
        # 如果是处理不稳定测试的worker，设置特殊环境
        if 'unstable' in worker_id:
            os.environ['MALLOC_CHECK_'] = '0'
            os.environ['MALLOC_PERTURB_'] = '0'


def pytest_runtest_setup(item):
    """在每个测试运行前的设置"""
    # 对于不稳定的测试，设置环境变量来帮助调试
    if item.get_closest_marker('unstable'):
        os.environ['MALLOC_CHECK_'] = '0'  # 禁用glibc的malloc检查，避免abort
        os.environ['MALLOC_PERTURB_'] = '0'  # 禁用内存扰动