"""
Pytest configuration and shared test utilities.

This file is automatically discovered by pytest and provides shared
utilities for all tests. Tests should prefer importing helpers from
``acceptance_log_utils`` to avoid depending on pytest's internal module
loading order.
"""

from pathlib import Path
import sys

_TESTS_DIR = Path(__file__).parent
if str(_TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(_TESTS_DIR))

from acceptance_log_utils import (  # noqa: F401
    acceptance_statistics,
    compute_detailed_balance_ratio,
    count_by_move_and_species,
    count_by_species,
    filter_by_move,
    filter_by_species,
    match_insert_delete_pairs,
    read_jsonl,
)

__all__ = [
    "read_jsonl",
    "count_by_species",
    "count_by_move_and_species",
    "filter_by_move",
    "filter_by_species",
    "match_insert_delete_pairs",
    "compute_detailed_balance_ratio",
    "acceptance_statistics",
]
