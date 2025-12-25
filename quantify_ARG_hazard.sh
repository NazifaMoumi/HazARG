#!/bin/bash

# Define sample IDs only
# sample_ids=($(ls -d /projects/ciwars/ARG_hazard/output/*/ | xargs -n 1 basename))
sample_ids=(
  SRR17257131
  SRR17257165
  SRR17257167
)

# Base paths
base_context_path="/projects/ciwars/ARG_context/ARGContextProfiler/results/murine"
base_output_path="/projects/ciwars/ARG_hazard/output/murine"

# Submit each sample as a separate Slurm job
for sample in "${sample_ids[@]}"; do
  context_fasta="${base_context_path}/${sample}/whole_context_with_length_clustered_filtered.fasta"
  output_dir="${base_output_path}/${sample}"

  sbatch --job-name="hazard_${sample}" \
         --output="/projects/ciwars/ARG_hazard/slurm_logs/murine_${sample}_%j.out" \
         --error="/projects/ciwars/ARG_hazard/slurm_logs/murine_${sample}_%j.err" \
         --time=80:00:00 \
         --partition=normal_q \
         --mem=100G \
         --account=metagen \
         --wrap="python -u /projects/ciwars/ARG_hazard/hazard_score_calculator.py --contexts_fasta '$context_fasta' --output_dir '$output_dir'"
done
