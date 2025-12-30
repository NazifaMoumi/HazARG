from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set

from .gene_id import GeneIdParser


def parse_mmseqs_tsv_unique_contexts(tsv_path: Path, parser: GeneIdParser) -> Dict[str, int]:
    """
    Count number of UNIQUE query context IDs aligning per gene key.

    The TSV must have at least:
      col0=query_id, col1=target_id
    """
    if not tsv_path.exists():
        return {}

    gene_to_contexts: Dict[str, Set[str]] = defaultdict(set)
    with tsv_path.open("r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row:
                continue
            query_id = row[0].strip()
            gene = parser.gene_key(query_id)
            gene_to_contexts[gene].add(query_id)

    return {g: len(s) for g, s in gene_to_contexts.items()}


def filter_self_hits_inplace(tsv_path: Path) -> None:
    """
    Optional: remove rows where query and target appear to be self-hits.
    This is a conservative heuristic that won't crash if formats differ.
    """
    if not tsv_path.exists():
        raise FileNotFoundError(tsv_path)

    kept = []
    with tsv_path.open("r", newline="") as infile:
        reader = csv.reader(infile, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            q, t = row[0], row[1]
            # heuristic: drop if the exact query id appears in target string
            if q and q in t:
                continue
            kept.append(row)

    with tsv_path.open("w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerows(kept)
