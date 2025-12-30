from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

from .config import MmseqsConfig, PipelineConfig, PresenceConfig, ThresholdConfig
from .pipeline import run_pipeline


def _parse_db_kv(pairs: list[str]) -> Dict[str, Path]:
    dbs: Dict[str, Path] = {}
    for p in pairs:
        if "=" not in p:
            raise ValueError(f"Invalid --db '{p}'. Use LABEL=/path/to/db.fasta")
        k, v = p.split("=", 1)
        dbs[k.strip()] = Path(v).expanduser().resolve()
    return dbs


def _parse_weights(pairs: list[str]) -> Dict[str, float]:
    w: Dict[str, float] = {}
    for p in pairs:
        if ":" not in p:
            raise ValueError(f"Invalid --weight '{p}'. Use LABEL:FLOAT")
        k, v = p.split(":", 1)
        w[k.strip()] = float(v)
    return w


def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="HazARG hazard index pipeline")
    ap.add_argument("--contexts-fasta", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)

    ap.add_argument("--query-gene", default="", help="If set, score only this gene key.")

    ap.add_argument("--lambda-low", type=float, default=10.0)
    ap.add_argument("--lambda-high", type=float, default=50.0)

    ap.add_argument("--cap-percentile", type=float, default=80.0)
    ap.add_argument("--low-pct", type=float, default=5.0)
    ap.add_argument("--high-pct", type=float, default=95.0)

    ap.add_argument("--mmseqs2", default="mmseqs")
    ap.add_argument("--mmseqs-total-threads", type=int, default=None)
    ap.add_argument("--mmseqs-max-jobs", type=int, default=None)

    ap.add_argument(
        "--db",
        action="append",
        default=[],
        help="Database mapping: LABEL=/path/to/db.fasta (repeatable).",
    )
    ap.add_argument(
        "--weight",
        action="append",
        default=[],
        help="Weight mapping: LABEL:FLOAT (repeatable). Defaults to 1.0 if absent.",
    )

    ap.add_argument("--verbose", action="store_true")
    return ap


def main() -> None:
    args = build_argparser().parse_args()

    databases = _parse_db_kv(args.db) if args.db else {}
    weights = _parse_weights(args.weight) if args.weight else {}

    cfg = PipelineConfig(
        contexts_fasta=args.contexts_fasta.expanduser().resolve(),
        output_dir=args.output_dir.expanduser().resolve(),
        query_gene=args.query_gene,
        thresholds=ThresholdConfig(
            cap_percentile=args.cap_percentile,
            low_pct=args.low_pct,
            high_pct=args.high_pct,
        ),
        presence=PresenceConfig(
            lambda_low=args.lambda_low,
            lambda_high=args.lambda_high,
        ),
        weights=weights if weights else PipelineConfig(contexts_fasta=Path("."), output_dir=Path(".")).weights,
        databases=databases,
        mmseqs=MmseqsConfig(
            mmseqs2_path=args.mmseqs2,
            total_threads=args.mmseqs_total_threads,
            max_jobs=args.mmseqs_max_jobs,
        ),
    )

    run_pipeline(cfg, verbose=args.verbose)
