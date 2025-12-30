from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class GeneIdParser:
    """
    Normalize headers/query IDs into a *gene key* used across the pipeline.

    Supports two formats commonly seen in your code:
    - pipe-delimited headers where you want: part1|part2|part3|part2
    - underscore-delimited headers like aadA_0012 -> aadA
    """

    def gene_key(self, raw_id: str) -> str:
        raw_id = raw_id.strip()
        if "|" in raw_id:
            parts = raw_id.split("|")
            if len(parts) >= 3:
                return "|".join([parts[0], parts[1], parts[2], parts[1]])
            # fallback: keep entire string if too short
            return raw_id
        # fallback underscore rule
        return raw_id.split("_", 1)[0]
