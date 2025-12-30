from __future__ import annotations

import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Tuple

from .config import MmseqsConfig


def _run(cmd: list[str]) -> Tuple[str, str]:
    proc = subprocess.run(
        cmd,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return proc.stdout, proc.stderr


def run_mmseqs_easy_search(
    *,
    contexts_fasta: Path,
    databases: Dict[str, Path],
    out_dir: Path,
    cfg: MmseqsConfig,
) -> Dict[str, Path]:
    """
    Run MMseqs2 easy-search against multiple target DBs in parallel.

    Returns:
        dict: db_label -> output_tsv_path
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    max_jobs = cfg.max_jobs or min(len(databases), (os.cpu_count() or 1))
    total_threads = cfg.total_threads or (os.cpu_count() or 1)
    threads_per_job = max(1, total_threads // max_jobs)

    commands: Dict[str, Tuple[list[str], Path]] = {}
    for label, db_path in databases.items():
        tmp_dir = out_dir / f"tmp_{label}"
        tmp_dir.mkdir(parents=True, exist_ok=True)

        out_tsv = out_dir / f"{label}_alignment.tsv"
        cmd = [
            cfg.mmseqs2_path, "easy-search",
            str(contexts_fasta), str(db_path), str(out_tsv), str(tmp_dir),
            "--search-type", str(cfg.search_type),
            "--format-output", cfg.format_output,
            "--threads", str(threads_per_job),
        ]
        commands[label] = (cmd, out_tsv)

    results: Dict[str, Path] = {}
    with ThreadPoolExecutor(max_workers=max_jobs) as ex:
        futs = {ex.submit(_run, cmd): label for label, (cmd, _) in commands.items()}
        for fut in as_completed(futs):
            label = futs[fut]
            try:
                fut.result()
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"MMseqs2 failed for {label}:\n{e.stderr}") from e
            results[label] = commands[label][1]

    return results
