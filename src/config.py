from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional


@dataclass(frozen=True)
class ThresholdConfig:
    cap_percentile: float = 80.0
    low_pct: float = 5.0
    high_pct: float = 95.0


@dataclass(frozen=True)
class PresenceConfig:
    lambda_low: float = 10.0
    lambda_high: float = 50.0


@dataclass(frozen=True)
class MmseqsConfig:
    mmseqs2_path: str = "mmseqs"
    search_type: str = "3"
    format_output: str = "query,target,evalue,qcov,tcov"
    max_jobs: Optional[int] = None          # default: number of DBs
    total_threads: Optional[int] = None     # default: os.cpu_count()


@dataclass(frozen=True)
class PipelineConfig:
    contexts_fasta: Path
    output_dir: Path

    # scoring knobs
    thresholds: ThresholdConfig = field(default_factory=ThresholdConfig)
    presence: PresenceConfig = field(default_factory=PresenceConfig)

    # weights for presence entities
    weights: Dict[str, float] = field(default_factory=lambda: {
        "ARGs_deepARG": 1.0,
        "MGEs_mobileOG": 1.0,
        "MRGs_BacMet": 1.0,
        "VFGs_VFDB": 1.0,
        "ESKAPEE": 1.0,
    })

    # DB label -> fasta/mmseqs db path (for easy-search targets)
    databases: Dict[str, Path] = field(default_factory=dict)

    # mmseqs runner config
    mmseqs: MmseqsConfig = field(default_factory=MmseqsConfig)

    # optional single-gene mode (otherwise scores all genes)
    query_gene: str = ""
