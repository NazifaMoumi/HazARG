# %%

import sys
import statistics
import matplotlib.pyplot as plt
import numpy as np
from typing import Dict, List, Tuple
import argparse
import os
import subprocess
from collections import defaultdict
import csv
import concurrent.futures
from Bio import SeqIO

def parse_arguments():
    """
    Parse command-line arguments and return an argparse.Namespace object.
    """
    parser = argparse.ArgumentParser(
        description="Hazard Scoring Pipeline for ARGs in Metagenomic Sample"
    )
    parser.add_argument(
        "--query_gene",
        # required=True,
        default='', # all genes from the context.fasta file
        help="Name (or ID) of the query ARG for which we want hazard scoring."
    )
    parser.add_argument(
        "--sample",
        # required=True,
        default = '',
        help="Identifier or path to the metagenomic sample for which we generate contexts."
    )
    parser.add_argument(
        "--argcontextprofiler_path",
        default="ARGContextProfiler",
        help="Path to the ARGContextProfiler tool."
    )
    parser.add_argument(
        "--mmseqs2_path",
        default="mmseqs",
        help="Path to the mmseqs2 binary."
    )
    parser.add_argument(
        "--output_dir",
        default="/projects/ciwars/ARG_hazard/output/test/",
        help="Directory to store intermediate and final results."
    )
    parser.add_argument(
        "--lambda_low",
        type=float,
        default=10.0,
        help="Lower percentage threshold (λ_low) for hazardous entities."
    )
    parser.add_argument(
        "--lambda_high",
        type=float,
        default=50.0,
        help="Upper percentage threshold (λ_high) for hazardous entities."
    )
    parser.add_argument(
        "--weights",
        nargs="+",
        default=["ARGs:1.0", "MGEs:1.0", "MRGs:1.0", "VFGs:1.0", "Pathogens:1.0"],
        help=(
            "List of weight assignments for each entity type in 'Entity:Weight' format. "
            "e.g. '--weights ARGs:2.0 MGEs:1.5 Pathogens:2.0'"
        )
    )
    parser.add_argument(
        "--cap_percentile",
        type=float,
        default=80.0,
        help="Percentile for capping outlier context counts (e.g., 80.0)."
    )
    parser.add_argument(
        "--low_pct",
        type=float,
        default=5.0,
        help="Percentile for capping outlier context counts (e.g., 80.0)."
    )
    parser.add_argument(
        "--high_pct",
        type=float,
        default=95.0,
        help="Percentile for capping outlier context counts (e.g., 80.0)."
    )
    parser.add_argument(
        "--contexts_fasta",
        default="/projects/ciwars/ARG_context/ARG-Context/BB_2_2021_01_08_INF_S61_analysis/result_contexts_cov_0.4/whole_context_with_length_clustered_filtered.fasta",
        help="Path to the FASTA file containing all genomic contexts (across multiple genes)."
    )
    return parser.parse_args()

def parse_alignment_tsv_per_gene(tsv_file: str) -> Dict[str, int]:
    """
    Parse an mmseqs2 TSV using the csv reader, counting how many UNIQUE contexts
    align for each gene. We assume:
       - The first column is the query_id (e.g. 'aadA_XXXX')
       - We parse the gene name by splitting on the first underscore ('aadA')
       - Each context (unique query_id) should be counted once per gene

    Return {gene -> number_of_unique_matched_contexts_in_this_DB}
    """

    if not os.path.exists(tsv_file):
        return {}

    gene_to_contexts = defaultdict(set)

    with open(tsv_file, "r", newline="") as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row:
                continue
            
            query_id = row[0].strip()
            parts = query_id.split("|")
            if len(parts) >= 3:
                # new name: first three parts + repeat of part 2
                gene = "|".join([parts[0], parts[1], parts[2], parts[1]])
            # Extract the gene name from the part before the first underscore
            # gene = query_id.split("_", 1)[0]
            gene_to_contexts[gene].add(query_id)

    gene_match_counts = {
        gene: len(context_ids)
        for gene, context_ids in gene_to_contexts.items()
    }

    return gene_match_counts

# def run_mmseqs2_align(contexts_fasta: str, databases: Dict[str, str], mmseqs2_path: str, output_dir: str) -> Dict[str, str]:
#     """
#     Run MMseqs2 alignments of the contexts FASTA against each reference database.

#     Args:
#         contexts_fasta (str): Path to the FASTA file containing genomic contexts.
#         databases (Dict[str, str]): Mapping of { 'NameOfDB': 'PathToMMseqsDB', ... }
#         mmseqs2_path (str): Path to mmseqs2 binary.
#         output_dir (str): Directory for storing alignment results.

#     Returns:
#         results_paths (Dict[str, str]): Mapping of DB name -> path to output TSV.
#     """
#     os.makedirs(output_dir, exist_ok=True)

#     results_paths = {}

#     for db_label, db_path in databases.items():
#         out_label = os.path.join(output_dir, f"{db_label}_alignment")
#         tsv_path = f"{out_label}.tsv"

#         # Example command: mmseqs easy-search <query> <targetDB> <result> <tmpDir> ...
#         command = [
#             mmseqs2_path, "easy-search",
#             contexts_fasta, db_path, tsv_path, os.path.join(output_dir, "tmp"), "--search-type","3",
#             "--format-output", "query,target,evalue,qcov,tcov"
#         ]
#         print(f"Running: {' '.join(command)}")
#         # commented out for testing 
#         # subprocess.run(command, check=True)

#         # The main results should be in out_label + ".tsv" if using easy-search
#         results_paths[db_label] = tsv_path

#     return results_paths

def run_mmseqs2_align(contexts_fasta: str,
                      databases: Dict[str, str],
                      mmseqs2_path: str,
                      output_dir: str) -> Dict[str, str]:
    """
    Run MMseqs2 alignments of the contexts FASTA against each reference database in parallel.
    """

    os.makedirs(output_dir, exist_ok=True)
    results_paths: Dict[str, str] = {}

    max_workers = min(len(databases), os.cpu_count() or 1)
    threads_per_job = max(1, (os.cpu_count() or 1) // max_workers)

    # Prepare commands
    commands = {}
    for db_label, db_path in databases.items():
        out_label = os.path.join(output_dir, f"{db_label}_alignment")
        tsv_path  = f"{out_label}.tsv"
        tmp_dir = os.path.join(output_dir, f"tmp_{db_label}")
        os.makedirs(tmp_dir, exist_ok=True)

        cmd = [
            mmseqs2_path, "easy-search",
            contexts_fasta, db_path, tsv_path, tmp_dir,
            "--search-type", "3",
            "--format-output", "query,target,evalue,qcov,tcov",
            "--threads", str(threads_per_job)
        ]
        commands[db_label] = (cmd, tsv_path)

    def _run_and_capture(db_label, cmd):
        try:
            # capture both stdout and stderr
            subprocess.run(cmd, check=True,
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           text=True)
        except subprocess.CalledProcessError as e:
            # print the stderr for debugging
            print(f"\n*** mmseqs2 failed for {db_label} ***\n{e.stderr}\n", flush=True)
            raise

    # Execute in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_label = {
            executor.submit(_run_and_capture, db_label, cmd): db_label
            for db_label, (cmd, _) in commands.items()
        }
        for fut in concurrent.futures.as_completed(future_to_label):
            label = future_to_label[fut]
            # this will re-raise CalledProcessError after printing stderr
            fut.result()
            results_paths[label] = commands[label][1]

    return results_paths

def calculate_presence_score_for_gene(
        gene: str,
        gene_match_counts: Dict[str, int],
        all_gene_counts: Dict[str, int],
        lambda_low: float,
        lambda_high: float
    ) -> int:
    """
    For a single gene:
      presence_pct = (# matched contexts for that gene / total contexts for that gene) * 100

    Then normalize that percentage to [0,1] by mapping:
      presence_pct <= lambda_low  -> 0.0
      presence_pct >= lambda_high -> 1.0
      linear interpolation in between.

    Returns:
      presence_score ∈ [0.0, 1.0]
    """
    total_contexts = all_gene_counts.get(gene, 0)
    if total_contexts == 0:
        return 0.0

    matched = gene_match_counts.get(gene, 0)
    presence_pct = (matched / total_contexts) * 100.0

    # continuous mapping:
    #   pct = lambda_low   -> 0.0
    #   pct = lambda_high  -> 1.0
    # beyond those bounds it’s clipped
    score = (presence_pct - lambda_low) / (lambda_high - lambda_low)
    return float(np.clip(score, 0.0, 1.0))

def compute_combined_score(entity_scores: Dict[str, int],
                           weights: Dict[str, float],
                           arg_hazard_score: int) -> float:
    """
    Sum up (entity_score * weight) + arg_hazard_score.
    """
    total = 0.0
    for entity_type, score_val in entity_scores.items():
        w = weights.get(entity_type, 1.0)
        total += score_val * w
    total += arg_hazard_score
    return total

def compute_thresholds_capped(
    gene_counts: Dict[str, int],
    cap_percentile: float = 80.0,
    low_pct: float = 5.0,
    high_pct: float = 95.0
) -> Tuple[float, float]:
    """
    1) Collect counts = list(gene_counts.values()).
    2) Winsorize/cap at the 'cap_percentile' (e.g. 80th) to reduce impact of extreme outliers.
    3) Compute tau_low as the low_pct-th and tau_high as the high_pct-th percentiles of the capped data.

    This gives you robust lower/upper bounds for scaling your context counts into [0,1].
    """
    if not gene_counts:
        raise ValueError("No gene counts found.")

    counts = np.array(list(gene_counts.values()), dtype=float)
    cap_val = np.percentile(counts, cap_percentile)
    capped = np.minimum(counts, cap_val)

    tau_low  = float(np.percentile(capped, low_pct))
    tau_high = float(np.percentile(capped, high_pct))

    # if everything is constant, avoid zero‐division later
    if tau_high <= tau_low:
        # fall back to min/max of capped
        tau_low, tau_high = float(capped.min()), float(capped.max())

    return tau_low, tau_high

def classify_arg_hazard(
    context_count: float,
    tau_low: float,
    tau_high: float
) -> float:
    """
    Map a gene's raw context_count into a continuous 0–1 hazard score:
      0 if count <= tau_low,
      1 if count >= tau_high,
      linear in between.
    """
    if context_count <= tau_low:
        return 0.0
    if context_count >= tau_high:
        return 1.0
    # linear interpolate
    return (context_count - tau_low) / (tau_high - tau_low)

def parse_context_counts(fasta_file: str) -> Dict[str, int]:
    """
    Returns a dict {gene -> number_of_contexts} for each gene in contexts.fasta.
    A header like '>aadA_0012' means gene='aadA'.
    """
    if not os.path.exists(fasta_file):
        raise FileNotFoundError(f"FASTA file {fasta_file} does not exist.")
    gene_counts = {}
    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:].strip()
                parts = header.split("|")
                if len(parts) >= 3:
                    # new name: first three parts + repeat of part 2
                    gene = "|".join([parts[0], parts[1], parts[2], parts[1]])
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
    return gene_counts

def filter_ARG_alignment_file(input_tsv):
    """
    Reads an mmseqs2 TSV using csv.reader (tab-delimited).
    Skips rows where the query gene name (the part before the first underscore in column 0)
    is found in the target field (column 1).
    """
    # output_tsv = '/projects/ciwars/ARG_hazard/output/mmseqs2_results/ARGs_deepARG_alignment_2.tsv'
    if not os.path.exists(input_tsv):
        raise FileNotFoundError(f"{input_tsv} not found.")

    filtered_rows = []
    
    with open(input_tsv, "r", newline="") as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            if len(row) < 2:
                continue
            query, target = row[0], row[1]
            query_gene_name = query.split("_", 1)[0]
            if query_gene_name in target:
                continue
            filtered_rows.append(row)

    # Overwrite the original file
    with open(input_tsv, "w", newline="") as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerows(filtered_rows)

def fix_gene_name(old_gene, contexts_fasta):
    """
    Given a gene identifier (partial) and a FASTA file, this function finds the best matching
    header in the FASTA file and constructs the full gene name as:
    part1|part2|part3|part2

    Parameters:
    - old_gene (str): Partial gene name like "FEATURES|bacterial"
    - contexts_fasta (str): Path to the FASTA file

    Returns:
    - new_gene_name (str): Updated gene name based on matched header
    """
    for record in SeqIO.parse(contexts_fasta, "fasta"):
        header = record.description
        if old_gene in header:
            parts = header.split("|")
            if len(parts) >= 3:
                # new name: first three parts + repeat of part 2
                new_name = "|".join([parts[0], parts[1], parts[2], parts[1]])
                return new_name
    return old_gene

def main():

    args = parse_arguments()

    # Convert the weights argument into a usable dict
    # e.g. ["ARGs:2.0", "MGEs:1.5"] -> {"ARGs":2.0, "MGEs":1.5}
    weight_dict = {}
    for w in args.weights:
        try:
            entity, val = w.split(":")
            weight_dict[entity] = float(val)
        except ValueError:
            print(f"Warning: Could not parse weight '{w}', skipping.")
            continue

    # Step 1: Run ARGContextProfiler to generate contexts FASTA
    # contexts_fasta = run_argcontextprofiler(
    #     query_gene=args.query_gene,
    #     sample=args.sample,
    #     output_dir=args.output_dir,
    #     profiler_path=args.argcontextprofiler_path
    # )

    print(f"\n[STEP] Parsing FASTA for context counts: {args.contexts_fasta}")
    all_gene_counts = parse_context_counts(args.contexts_fasta)  # { gene: count_of_contexts }
    if not all_gene_counts:
        print("No gene contexts found. Exiting.")
        sys.exit(1)

    print(f"Found {len(all_gene_counts)} distinct genes in {args.contexts_fasta}")

    # B) Compute global thresholds (tau_low, tau_high) from the distribution
    print("[STEP] Computing (tau_low, tau_high) from context distribution (with capping).")
    tau_low, tau_high = compute_thresholds_capped(all_gene_counts, args.cap_percentile, args.low_pct, args.high_pct)
    print(f"   τ_low = {tau_low:.2f}, τ_high = {tau_high:.2f}")

    # query_gene_count = all_gene_counts.get(args.query_gene, 0)
    # query_arg_hazard_score = classify_arg_hazard(query_gene_count, tau_low, tau_high)

    # # Show that classification
    # print(f"Query gene: {args.query_gene}, context count = {query_gene_count}")
    # print(f"Query gene hazard score => {query_arg_hazard_score}")

    databases = {
        "ARGs_deepARG": "/projects/ciwars/databases/deeparg/dataset.fasta",
        "MGEs_mobileOG": "/projects/ciwars/databases/mobileOG-db_beatrix-1.6.All.faa",
        "MRGs_BacMet": "/projects/ciwars/databases/BacMet2_EXP_database.fasta",
        "VFGs_VFDB": "/projects/ciwars/databases/VFDB_setA_pro.fas",
        # "Pathogens_GTDB": "/path/to/GTDB",
        "ESKAPEE": "/projects/ciwars/ARG_hazard/merged_ESKAPEE.fasta"
    }
    
    os.makedirs(args.output_dir, exist_ok=True)

    # Step 3: Run MMseqs2 searches
    alignment_results = run_mmseqs2_align(
        contexts_fasta=args.contexts_fasta,
        databases=databases,
        mmseqs2_path=args.mmseqs2_path,
        output_dir=os.path.join(args.output_dir, "mmseqs2_results")
    )
    alignment_results = {
        db_label: os.path.join(args.output_dir, "mmseqs2_results", f"{db_label}_alignment.tsv")
        for db_label in databases
    }
    filter_ARG_alignment_file(alignment_results['ARGs_deepARG'])

    alignment_matches = {}
    for db_label, tsv_path in alignment_results.items():
        gene_match_counts = parse_alignment_tsv_per_gene(tsv_path)
        alignment_matches[db_label] = gene_match_counts

    results = []
    for gene, context_count in all_gene_counts.items():

        arg_hazard_score = classify_arg_hazard(context_count, tau_low, tau_high)
        entity_scores = {}
        for db_label in alignment_matches.keys():
            entity_score = calculate_presence_score_for_gene(
                gene=gene,
                gene_match_counts=alignment_matches[db_label],
                all_gene_counts=all_gene_counts,
                lambda_low=args.lambda_low,
                lambda_high=args.lambda_high
            )
            entity_scores[db_label] = entity_score
        final_score = compute_combined_score(entity_scores, weight_dict, arg_hazard_score)
        final_score = final_score / (sum(weight_dict.values()) + 1.0)

        # new_gene = fix_gene_name(gene, args.contexts_fasta)
        row = {
            "Gene": gene,
            "ContextCount": context_count,
            "HazardScore": arg_hazard_score,
            **{f"{db_label}_PresenceScore": entity_scores[db_label] for db_label in entity_scores},
            "CombinedScore": final_score
        }
        results.append(row)

    # Sort results by final hazard or gene name, optionally
    results.sort(key=lambda x: x["CombinedScore"], reverse=True)

    # F) Print results
    print("\n========== HAZARD SCORING RESULTS (ALL GENES) ==========")
    for row in results:
        gene = row["Gene"]
        hazard = row["HazardScore"]
        combined = row["CombinedScore"]
        print(f"Gene: {gene}  |  CtxCount={row['ContextCount']}  "
              f"Hazard={hazard}  Combined={combined:.2f}")

    # G) Optionally, write to CSV
    out_csv = os.path.join(args.output_dir, "hazard_scores_all_genes.csv")
    fieldnames = list(results[0].keys())
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)
    print(f"[INFO] Saved detailed results to {out_csv}")

if __name__ == "__main__":
    main()

# %%