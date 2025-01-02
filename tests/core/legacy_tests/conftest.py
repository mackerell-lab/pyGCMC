import pytest
from pathlib import Path

def pytest_ignore_collect(collection_path: Path, config):
    """忽略legacy_tests目录下的所有测试文件"""
    return "legacy_tests" in str(collection_path) 