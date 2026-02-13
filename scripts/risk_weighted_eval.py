#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

CALLERS = ["gatk", "deepvariant", "strelka2", "freebayes"]
DATASETS = ["raw", "wes_standard"]
MODELS = ["consequence", "alphamissense", "alphagenome", "varsage", "max_all"]

CONSEQUENCE_WEIGHTS = {
    "transcript_ablation": 1.0,
    "splice_acceptor_variant": 1.0,
    "splice_donor_variant": 1.0,
    "stop_gained": 0.95,
    "frameshift_variant": 0.95,
    "stop_lost": 0.9,
    "start_lost": 0.9,
    "inframe_insertion": 0.7,
    "inframe_deletion": 0.7,
    "missense_variant": 0.6,
    "protein_altering_variant": 0.6,
    "splice_region_variant": 0.5,
    "synonymous_variant": 0.1,
    "intron_variant": 0.1,
    "upstream_gene_variant": 0.05,
    "downstream_gene_variant": 0.05,
    "intergenic_variant": 0.01,
}


Variant = Tuple[str, str, str, str, str]
Key = Tuple[str, str, str, str]


def parse_args():
    p = argparse.ArgumentParser(description="Risk-weighted FP/FN evaluation across models and WES filtering")
    p.add_argument("--variants-dir", required=True)
    p.add_argument("--prefix", required=True)
    p.add_argument("--errors-dir", required=True)
    p.add_argument("--output-summary", required=True)
    p.add_argument("--output-details", required=True)
    p.add_argument("--alpha-missense", default="")
    p.add_argument("--alpha-genome", default="")
    p.add_argument("--varsage", default="")
    return p.parse_args()


def load_score_table(path: str) -> Dict[Key, float]:
    scores: Dict[Key, float] = {}
    if not path:
        return scores
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return scores

    with p.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row.get("CHROM", ""), row.get("POS", ""), row.get("REF", ""), row.get("ALT", ""))
            try:
                scores[key] = float(row.get("SCORE", ""))
            except ValueError:
                continue
    return scores


def parse_info_for_consequence(info: str) -> float:
    consequence_weight = 0.2
    for field in info.split(";"):
        if field.startswith("ANN=") or field.startswith("CSQ="):
            ann = field.split("=", 1)[1]
            effects = []
            for rec in ann.split(","):
                parts = rec.split("|")
                if len(parts) > 1:
                    effects.extend(parts[1].split("&"))
            if effects:
                consequence_weight = max(CONSEQUENCE_WEIGHTS.get(e, 0.2) for e in effects)
            break
    return consequence_weight


def iter_vcf(path: Path) -> Iterable[Variant]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, flt, info = cols[:8]
            for a in alt.split(","):
                yield chrom, pos, ref, a, info


def load_vcf_variants(path: Path) -> List[Variant]:
    if not path.exists():
        return []
    return list(iter_vcf(path))


def weight_for_model(
    key: Key,
    info: str,
    model: str,
    am_scores: Dict[Key, float],
    ag_scores: Dict[Key, float],
    vs_scores: Dict[Key, float],
) -> float:
    consequence = parse_info_for_consequence(info)
    am = am_scores.get(key, 0.0)
    ag = ag_scores.get(key, 0.0)
    vs = vs_scores.get(key, 0.0)

    if model == "consequence":
        val = consequence
    elif model == "alphamissense":
        val = am if am > 0 else consequence
    elif model == "alphagenome":
        val = ag if ag > 0 else consequence
    elif model == "varsage":
        val = vs if vs > 0 else consequence
    else:  # max_all
        val = max(consequence, am, ag, vs)

    return min(max(val, 0.01), 1.0)


def summarize_set(
    variants: List[Variant],
    model: str,
    am_scores: Dict[Key, float],
    ag_scores: Dict[Key, float],
    vs_scores: Dict[Key, float],
):
    count = 0
    total_weight = 0.0
    high_impact = 0
    for chrom, pos, ref, alt, info in variants:
        count += 1
        key = (chrom, pos, ref, alt)
        w = weight_for_model(key, info, model, am_scores, ag_scores, vs_scores)
        total_weight += w
        if w >= 0.7:
            high_impact += 1
    return count, total_weight, high_impact


def main():
    args = parse_args()

    am_scores = load_score_table(args.alpha_missense)
    ag_scores = load_score_table(args.alpha_genome)
    vs_scores = load_score_table(args.varsage)

    details_rows = []
    summary_rows = []

    for caller in CALLERS:
        caller_dir = Path(args.variants_dir) / caller

        pass_raw = caller_dir / f"{args.prefix}_{caller}_pass.norm.vcf.gz"
        pass_wes = caller_dir / f"{args.prefix}_{caller}_wes_standard.pass.norm.vcf.gz"

        dataset_to_pass = {
            "raw": pass_raw,
            "wes_standard": pass_wes,
        }

        for dataset in DATASETS:
            pass_path = dataset_to_pass[dataset]
            if not pass_path.exists():
                continue

            fp_vcf = Path(args.errors_dir) / caller / dataset / f"{args.prefix}_{caller}_{dataset}_fp.vcf.gz"
            fn_vcf = Path(args.errors_dir) / caller / dataset / f"{args.prefix}_{caller}_{dataset}_fn.vcf.gz"

            pass_variants = load_vcf_variants(pass_path)
            fp_variants = load_vcf_variants(fp_vcf)
            fn_variants = load_vcf_variants(fn_vcf)

            for model in MODELS:
                pass_count, pass_weight, _ = summarize_set(pass_variants, model, am_scores, ag_scores, vs_scores)
                fp_count, fp_weight, fp_hi = summarize_set(fp_variants, model, am_scores, ag_scores, vs_scores)
                fn_count, fn_weight, fn_hi = summarize_set(fn_variants, model, am_scores, ag_scores, vs_scores)

                tp_weight = max(pass_weight - fp_weight, 0.0)
                tp_count = max(pass_count - fp_count, 0)

                w_precision = tp_weight / (tp_weight + fp_weight) if (tp_weight + fp_weight) > 0 else 0.0
                w_recall = tp_weight / (tp_weight + fn_weight) if (tp_weight + fn_weight) > 0 else 0.0
                w_f1 = (2 * w_precision * w_recall / (w_precision + w_recall)) if (w_precision + w_recall) > 0 else 0.0

                details_rows.append({
                    "Caller": caller,
                    "Dataset": dataset,
                    "Model": model,
                    "FP_Count": fp_count,
                    "FN_Count": fn_count,
                    "FP_WeightedCost": f"{fp_weight:.6f}",
                    "FN_WeightedCost": f"{fn_weight:.6f}",
                    "FP_HighImpact_n": fp_hi,
                    "FN_HighImpact_n": fn_hi,
                })

                summary_rows.append({
                    "Caller": caller,
                    "Dataset": dataset,
                    "Model": model,
                    "TP_Count": tp_count,
                    "TP_Weighted": f"{tp_weight:.6f}",
                    "FP_Weighted": f"{fp_weight:.6f}",
                    "FN_Weighted": f"{fn_weight:.6f}",
                    "WeightedPrecision": f"{w_precision:.6f}",
                    "WeightedRecall": f"{w_recall:.6f}",
                    "WeightedF1": f"{w_f1:.6f}",
                })

    Path(args.output_details).parent.mkdir(parents=True, exist_ok=True)

    with Path(args.output_details).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Caller", "Dataset", "Model", "FP_Count", "FN_Count",
                "FP_WeightedCost", "FN_WeightedCost", "FP_HighImpact_n", "FN_HighImpact_n",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(details_rows)

    with Path(args.output_summary).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Caller", "Dataset", "Model", "TP_Count", "TP_Weighted", "FP_Weighted",
                "FN_Weighted", "WeightedPrecision", "WeightedRecall", "WeightedF1",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)


if __name__ == "__main__":
    main()
