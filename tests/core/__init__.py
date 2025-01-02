import os
import sys
import shutil
import pytest
from pathlib import Path

# 设置Python不生成字节码文件
sys.dont_write_bytecode = True
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'

def clean_pycache():
    """清理所有__pycache__目录"""
    test_root = os.path.dirname(__file__)
    cleaned = False
    for root, dirs, files in os.walk(test_root):
        for dir_name in dirs:
            if dir_name == "__pycache__":
                pycache_path = os.path.join(root, dir_name)
                try:
                    shutil.rmtree(pycache_path)
                    print(f"Cleaned: {pycache_path}")
                    cleaned = True
                except Exception as e:
                    print(f"Failed to clean {pycache_path}: {e}")
    if not cleaned:
        print("No __pycache__ directories found to clean")

@pytest.fixture(scope="session", autouse=True)
def cleanup_after_tests(request):
    """测试会话结束后自动清理__pycache__"""
    def cleanup():
        print("\nCleaning up __pycache__ directories...")
        clean_pycache()
    request.addfinalizer(cleanup)

def pytest_sessionfinish(session, exitstatus):
    """在pytest会话结束时清理__pycache__"""
    print("\nTest session finished, cleaning up...")
    clean_pycache()

def pytest_ignore_collect(collection_path: Path, config):
    """全局忽略legacy_tests目录"""
    return "legacy_tests" in str(collection_path) 