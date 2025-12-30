from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np


@dataclass(frozen=True)
class Thresholds:
    tau_low: float
    tau_high: float


def compute_thresholds_capped(
    gene_counts: Dict[str, int],
    *,
    cap_percentile: float,
    low_pct: float,
    high_pct: float,
) -> Thresholds:
    """
    Winsorize counts at `cap_percentile`, then compute tau_low/tau_high as
    percentiles of capped values.
    """
    if not gene_counts:
        raise ValueError("No gene counts found.")

    counts = np.array(list(gene_counts.values()), dtype=float)
    cap_val = np.percentile(counts, cap_percentile)
    capped = np.minimum(counts, cap_val)

    tau_low = float(np.percentile(capped, low_pct))
    tau_high = float(np.percentile(capped, high_pct))
    if tau_high <= tau_low:
        tau_low, tau_high = float(capped.min()), float(capped.max())

    return Thresholds(tau_low=tau_low, tau_high=tau_high)


def hazard_index_from_context_count(context_count: float, thr: Thresholds) -> float:
    """Map a context count into [0,1] using (tau_low, tau_high)."""
    if context_count <= thr.tau_low:
        return 0.0
    if context_count >= thr.tau_high:
        return 1.0
    return (context_count - thr.tau_low) / (thr.tau_high - thr.tau_low)


def presence_score(
    *,
    gene: str,
    matched_contexts: Dict[str, int],
    total_contexts: Dict[str, int],
    lambda_low: float,
    lambda_high: float,
) -> float:
    """
    presence_pct = matched/total * 100
    then linearly map into [0,1] between lambda_low and lambda_high.
    """
    tot = total_contexts.get(gene, 0)
    if tot == 0:
        return 0.0
    matched = matched_contexts.get(gene, 0)
    pct = (matched / tot) * 100.0
    score = (pct - lambda_low) / (lambda_high - lambda_low)
    return float(np.clip(score, 0.0, 1.0))


def combined_index(
    *,
    hazard_index: float,
    entity_presence: Dict[str, float],
    weights: Dict[str, float],
) -> float:
    """
    Weighted mean of:
      hazard_index (weight 1.0)
      + sum(entity_presence[label] * weights[label])

    Normalized by (1 + sum(weights)).
    """
    total = hazard_index
    denom = 1.0
    for label, val in entity_presence.items():
        w = float(weights.get(label, 1.0))
        total += val * w
        denom += w
    return total / denom
