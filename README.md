# HazARG
## Hazard Index Pipeline

Compute per-gene hazard index from:
1) Context multiplicity (winsorized percentile thresholds)
2) Presence scores from MMseqs2 alignments vs multiple entity DBs
3) Weighted aggregation into a final combined index

## Install (editable)
pip install -e .

## Run
hazarg-score \
  --contexts-fasta /path/to/contexts.fasta \
  --output-dir /path/to/out \
  --mmseqs2 mmseqs \
  --db ARGs_deepARG=/db/deeparg.fasta \
  --db MGEs_mobileOG=/db/mobileOG.faa \
  --db MRGs_BacMet=/db/bacmet.fasta \
  --db VFGs_VFDB=/db/vfdb.faa \
  --db ESKAPEE=/db/eskapee.fasta

## Output
- out/hazard_scores_all_genes.csv
- out/logs/run.log
- out/mmseqs2_results/*.tsv
