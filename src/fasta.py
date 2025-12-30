from __future__ import annotations

from pathlib import Path
from typing import Dict

from Bio import SeqIO

from .gene_id import GeneIdParser


def count_contexts_by_gene(contexts_fasta: Path, parser: GeneIdParser) -> Dict[str, int]:
    """
    Count how many context sequences exist per gene.

    Returns:
        dict: gene_key -> count
    """
    if not contexts_fasta.exists():
        raise FileNotFoundError(f"FASTA not found: {contexts_fasta}")

    counts: Dict[str, int] = {}
    with contexts_fasta.open("r") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            gene = parser.gene_key(rec.description)
            counts[gene] = counts.get(gene, 0) + 1
    return counts
