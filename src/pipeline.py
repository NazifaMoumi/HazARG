from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Dict, List

from .config import PipelineConfig
from .fasta import count_contexts_by_gene
from .gene_id import GeneIdParser
from .io_utils import write_dict_rows_csv
from .logging_utils import setup_logging
from .mmseqs import run_mmseqs_easy_search
from .parsing import filter_self_hits_inplace, parse_mmseqs_tsv_unique_contexts
from .scoring import (
    combined_index,
    compute_thresholds_capped,
    hazard_index_from_context_count,
    presence_score,
)


def run_pipeline(cfg: PipelineConfig, *, verbose: bool = False) -> Path:
    logger = setup_logging(cfg.output_dir, verbose=verbose)
    parser = GeneIdParser()

    logger.info(f"Parsing contexts: {cfg.contexts_fasta}")
    gene_counts = count_contexts_by_gene(cfg.contexts_fasta, parser)
    if not gene_counts:
        raise RuntimeError("No gene contexts found in FASTA.")

    logger.info(f"Found {len(gene_counts)} genes with contexts.")

    thr = compute_thresholds_capped(
        gene_counts,
        cap_percentile=cfg.thresholds.cap_percentile,
        low_pct=cfg.thresholds.low_pct,
        high_pct=cfg.thresholds.high_pct,
    )
    logger.info(f"Thresholds: tau_low={thr.tau_low:.2f}, tau_high={thr.tau_high:.2f}")

    if not cfg.databases:
        raise ValueError("No databases provided. Use --db LABEL=/path/to/db.fasta")

    mm_out = cfg.output_dir / "mmseqs2_results"
    logger.info(f"Running MMseqs2 vs {len(cfg.databases)} DBs -> {mm_out}")
    tsvs = run_mmseqs_easy_search(
        contexts_fasta=cfg.contexts_fasta,
        databases=cfg.databases,
        out_dir=mm_out,
        cfg=cfg.mmseqs,
    )

    # Optional cleanup for a specific DB (kept as generic self-hit filter)
    if "ARGs_deepARG" in tsvs:
        logger.info("Filtering possible self-hits in ARGs_deepARG TSV")
        filter_self_hits_inplace(tsvs["ARGs_deepARG"])

    logger.info("Parsing MMseqs2 TSVs into per-gene unique context-hit counts")
    match_counts_by_db: Dict[str, Dict[str, int]] = {
        label: parse_mmseqs_tsv_unique_contexts(path, parser)
        for label, path in tsvs.items()
    }

    genes = [cfg.query_gene] if cfg.query_gene else list(gene_counts.keys())

    rows: List[Dict] = []
    for gene in genes:
        ctx = float(gene_counts.get(gene, 0))
        hz = hazard_index_from_context_count(ctx, thr)

        entity_presence: Dict[str, float] = {}
        for db_label, matched_counts in match_counts_by_db.items():
            entity_presence[db_label] = presence_score(
                gene=gene,
                matched_contexts=matched_counts,
                total_contexts=gene_counts,
                lambda_low=cfg.presence.lambda_low,
                lambda_high=cfg.presence.lambda_high,
            )

        comb = combined_index(hazard_index=hz, entity_presence=entity_presence, weights=cfg.weights)

        row = {
            "Gene": gene,
            "ContextCount": int(ctx),
            "HazardIndex": hz,
            **{f"{k}_PresenceScore": v for k, v in entity_presence.items()},
            "CombinedIndex": comb,
        }
        rows.append(row)

    rows.sort(key=lambda r: r["CombinedIndex"], reverse=True)

    out_csv = cfg.output_dir / "hazard_scores_all_genes.csv"
    write_dict_rows_csv(out_csv, rows)
    logger.info(f"Wrote: {out_csv}")
    return out_csv
