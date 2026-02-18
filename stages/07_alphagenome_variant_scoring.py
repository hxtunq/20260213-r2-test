#!/usr/bin/env python3
"""Score FN/FP variants from rtg vcfeval outputs with AlphaGenome on terminal Ubuntu."""

from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path
from typing import Iterable



DEFAULT_CALLERS = ["gatk", "deepvariant", "freebayes", "strelka2"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Load fn/fp.vcf.gz from rtg vcfeval outputs for multiple callers and "
            "score variants with AlphaGenome."
        )
    )
    parser.add_argument(
        "--rtg-root",
        type=Path,
        default=Path("/workspace/20260213-r2-test/results/benchmarks/rtg_vcfeval"),
        help="Root directory containing caller subdirs with fn.vcf.gz/fp.vcf.gz.",
    )
    parser.add_argument(
        "--callers",
        nargs="+",
        default=DEFAULT_CALLERS,
        help="Caller subdirectories to read.",
    )
    parser.add_argument(
        "--variant-classes",
        nargs="+",
        choices=["fn", "fp"],
        default=["fn", "fp"],
        help="Variant classes to score.",
    )
    parser.add_argument(
        "--organism",
        choices=["human", "mouse"],
        default="human",
        help="Organism for AlphaGenome scoring.",
    )
    parser.add_argument(
        "--sequence-length",
        choices=["16KB", "100KB", "500KB", "1MB"],
        default="1MB",
        help="Sequence context length around each variant.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("variant_scores_fn-fp_4callers.csv"),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--count-only",
        action="store_true",
        help="Only count and print variants found, do not call AlphaGenome API.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Optional max variants to score (0 means all).",
    )
    return parser.parse_args()


def load_vcf_as_variants(vcf_path: Path, caller: str, variant_class: str) -> "pd.DataFrame":
    import pandas as pd

    rows = []
    with gzip.open(vcf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom, pos, _id, ref, alt = fields[:5]
            for alt_allele in alt.split(","):
                rows.append(
                    {
                        "variant_id": f"{caller}_{variant_class}_{chrom}_{pos}_{ref}_{alt_allele}",
                        "CHROM": chrom,
                        "POS": int(pos),
                        "REF": ref,
                        "ALT": alt_allele,
                        "caller": caller,
                        "variant_class": variant_class,
                    }
                )
    return pd.DataFrame(rows)


def load_variants(rtg_root: Path, callers: Iterable[str], variant_classes: Iterable[str]) -> "pd.DataFrame":
    import pandas as pd

    frames = []
    for variant_class in variant_classes:
        for caller in callers:
            vcf_path = rtg_root / caller / f"{variant_class}.vcf.gz"
            if not vcf_path.exists():
                print(f"[WARN] Missing file: {vcf_path}")
                continue
            frame = load_vcf_as_variants(vcf_path, caller, variant_class)
            if frame.empty:
                print(f"[WARN] Empty variants: {vcf_path}")
                continue
            frames.append(frame)

    if not frames:
        raise ValueError("No variants were loaded from provided rtg-root/callers/classes.")

    vcf = pd.concat(frames, ignore_index=True)
    required_columns = ["variant_id", "CHROM", "POS", "REF", "ALT"]
    missing = [c for c in required_columns if c not in vcf.columns]
    if missing:
        raise ValueError(f"VCF-derived table missing required columns: {missing}")

    return vcf




def count_variants(rtg_root: Path, callers: Iterable[str], variant_classes: Iterable[str]) -> int:
    total = 0
    for variant_class in variant_classes:
        for caller in callers:
            vcf_path = rtg_root / caller / f"{variant_class}.vcf.gz"
            if not vcf_path.exists():
                print(f"[WARN] Missing file: {vcf_path}")
                continue
            with gzip.open(vcf_path, "rt") as fh:
                count = sum(1 for line in fh if line and not line.startswith("#"))
            print(f"{caller}	{variant_class}	{count}")
            total += count
    print(f"Total variants loaded: {total:,}")
    return total

def score_with_alphagenome(vcf: pd.DataFrame, organism: str, sequence_length: str) -> "pd.DataFrame":
    from alphagenome import colab_utils
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers
    from tqdm import tqdm

    dna_model = dna_client.create(colab_utils.get_api_key())

    organism_map = {
        "human": dna_client.Organism.HOMO_SAPIENS,
        "mouse": dna_client.Organism.MUS_MUSCULUS,
    }
    organism_enum = organism_map[organism]
    seq_len = dna_client.SUPPORTED_SEQUENCE_LENGTHS[f"SEQUENCE_LENGTH_{sequence_length}"]

    selected_scorers = list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values())

    unsupported = [
        scorer
        for scorer in selected_scorers
        if (
            organism_enum.value
            not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
        )
        or (
            scorer.requested_output == dna_client.OutputType.PROCAP
            and organism_enum == dna_client.Organism.MUS_MUSCULUS
        )
    ]
    for scorer in unsupported:
        selected_scorers.remove(scorer)

    results = []
    for _, row in tqdm(vcf.iterrows(), total=len(vcf)):
        variant = genome.Variant(
            chromosome=str(row.CHROM),
            position=int(row.POS),
            reference_bases=row.REF,
            alternate_bases=row.ALT,
            name=row.variant_id,
        )
        interval = variant.reference_interval.resize(seq_len)
        results.append(
            dna_model.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=selected_scorers,
                organism=organism_enum,
            )
        )

    df_scores = variant_scorers.tidy_scores(results)
    metadata = vcf[["variant_id", "caller", "variant_class"]].drop_duplicates()
    return df_scores.merge(metadata, on="variant_id", how="left")


def main() -> None:
    args = parse_args()

    if args.limit < 0:
        raise ValueError("--limit must be >= 0")

    if args.count_only:
        count_variants(args.rtg_root, args.callers, args.variant_classes)
        return

    vcf = load_variants(args.rtg_root, args.callers, args.variant_classes)

    print(vcf.groupby(["caller", "variant_class"]).size())
    print(f"Total variants loaded: {len(vcf):,}")

    if args.limit > 0:
        vcf = vcf.head(args.limit).copy()
        print(f"Applying --limit, variants kept for scoring: {len(vcf):,}")

    if "ALPHAGENOME_API_KEY" not in os.environ:
        print("[INFO] ALPHAGENOME_API_KEY is not set. You may be prompted by SDK flow.")

    df_scores = score_with_alphagenome(vcf, args.organism, args.sequence_length)
    df_scores.to_csv(args.output, index=False)
    print(f"Saved scores: {args.output}")


if __name__ == "__main__":
    main()
