"""
Shared utilities for parsing and analyzing GCMC acceptance logs.

These helpers are deliberately kept outside of any pytest ``conftest`` module
so they can be imported reliably from tests without depending on pytest's
module-loading order.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Tuple

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


def read_jsonl(filepath: Path) -> List[Dict[str, Any]]:
    """Read JSONL file and return list of records."""
    records: List[Dict[str, Any]] = []

    if not filepath.exists():
        return records

    with open(filepath, "r") as handle:
        for line in handle:
            data = line.strip()
            if not data:
                continue
            try:
                records.append(json.loads(data))
            except json.JSONDecodeError as exc:
                print(
                    f"Warning: Failed to parse line: {data[:50]}... Error: {exc}"
                )

    return records


def count_by_species(records: List[Dict[str, Any]]) -> Dict[str, int]:
    """Count total attempts per species."""
    counts: Dict[str, int] = defaultdict(int)
    for record in records:
        species = record.get("species", "unknown")
        counts[species] += 1
    return dict(counts)


def count_by_move_and_species(
    records: List[Dict[str, Any]]
) -> Dict[Tuple[str, str], int]:
    """Count attempts grouped by (move_type, species)."""
    counts: Dict[Tuple[str, str], int] = defaultdict(int)
    for record in records:
        move = record.get("move", "unknown")
        species = record.get("species", "unknown")
        counts[(move, species)] += 1
    return dict(counts)


def filter_by_move(
    records: List[Dict[str, Any]], move_type: str
) -> List[Dict[str, Any]]:
    """Filter acceptance records by move type."""
    return [record for record in records if record.get("move") == move_type]


def filter_by_species(
    records: List[Dict[str, Any]], species: str
) -> List[Dict[str, Any]]:
    """Filter acceptance records by species name."""
    return [record for record in records if record.get("species") == species]


def match_insert_delete_pairs(
    records: List[Dict[str, Any]],
    species: str,
    max_step_gap: int = 100,
) -> List[Tuple[Dict[str, Any], Dict[str, Any]]]:
    """
    Match insertion-deletion pairs for detailed balance verification.

    Pairs are matched by species and by ensuring their MC steps are within
    ``max_step_gap`` of one another.
    """
    insertions = [
        record
        for record in records
        if record.get("move") == "insertion" and record.get("species") == species
    ]
    deletions = [
        record
        for record in records
        if record.get("move") == "deletion" and record.get("species") == species
    ]

    pairs: List[Tuple[Dict[str, Any], Dict[str, Any]]] = []
    used_deletions = set()

    for ins in insertions:
        ins_step = ins.get("step", 0)
        best_del: Tuple[int, Dict[str, Any]] | None = None
        best_gap = max_step_gap + 1

        for idx, deletion in enumerate(deletions):
            if idx in used_deletions:
                continue

            gap = abs(deletion.get("step", 0) - ins_step)
            if gap < best_gap:
                best_gap = gap
                best_del = (idx, deletion)

        if best_del is not None:
            used_deletions.add(best_del[0])
            pairs.append((ins, best_del[1]))

    return pairs


def compute_detailed_balance_ratio(
    insertion: Dict[str, Any], deletion: Dict[str, Any]
) -> float:
    """Compute π(A)T(A→B) / π(B)T(B→A) from paired moves."""
    p_ins = insertion.get("pAcc", 0.0)
    p_del = deletion.get("pAcc", 0.0)

    q_ratio_ins = insertion.get("proposalRatio", 1.0)
    q_ratio_del = deletion.get("proposalRatio", 1.0)

    if p_del > 0:
        return p_ins / p_del * q_ratio_del / q_ratio_ins
    return float("inf") if p_ins > 0 else 1.0


def acceptance_statistics(
    records: List[Dict[str, Any]]
) -> Dict[str, Dict[str, float]]:
    """Compute acceptance rates grouped by move and species."""
    stats: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(
        lambda: defaultdict(lambda: {"attempts": 0, "accepted": 0})
    )

    for record in records:
        move = record.get("move", "unknown")
        species = record.get("species", "unknown")
        accepted = record.get("accepted", False)

        stats[move][species]["attempts"] += 1
        if accepted:
            stats[move][species]["accepted"] += 1

    result: Dict[str, Dict[str, float]] = {}
    for move, species_dict in stats.items():
        result[move] = {}
        for species, counts in species_dict.items():
            attempts = counts["attempts"]
            result[move][species] = (
                counts["accepted"] / attempts if attempts > 0 else 0.0
            )

    return result

