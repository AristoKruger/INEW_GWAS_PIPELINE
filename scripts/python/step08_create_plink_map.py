import argparse
from pathlib import Path

import pandas as pd


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create a PLINK .map file by combining the filtered SNP list with the probe metadata",
    )
    parser.add_argument(
        "--snp-list",
        default="outputs/qc/filtered_snp_list.txt",
        type=Path,
        help="Path to the text file containing SNP IDs (one per line)",
    )
    parser.add_argument(
        "--metadata",
        default="data/metadata/35k_probe_set_IWGSCv1_cleaned_tab.map",
        type=Path,
        help="Tab-delimited annotation file with columns 35K SNPId, IWGSC_v1_Chromosome, IWGSC_v1_Position",
    )
    parser.add_argument(
        "--output",
        default="outputs/plink/filtered_genotypes_strict.map",
        type=Path,
        help="Destination for the PLINK .map output",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    snp_list = pd.read_csv(args.snp_list, header=None, names=["SNP"])
    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)

    meta = meta[meta["IWGSC_v1_Chromosome"].notna()]
    meta["IWGSC_v1_Position"] = pd.to_numeric(meta["IWGSC_v1_Position"], errors="coerce")
    meta = meta.dropna(subset=["IWGSC_v1_Position"])
    meta["Chrom"] = meta["IWGSC_v1_Chromosome"].str.replace("chr", "", regex=False)
    meta["Position"] = meta["IWGSC_v1_Position"].astype(int)

    merged = snp_list.merge(meta, left_on="SNP", right_on="35K SNPId", how="left")
    missing = merged[merged["Chrom"].isna()]
    print(f"Missing SNPs: {missing.shape[0]}")

    plink_map = merged[["Chrom", "SNP", "IWGSC_v1_Position"]].copy()
    plink_map.insert(2, "cM", 0)
    plink_map.columns = [0, 1, 2, 3]
    plink_map = plink_map.dropna()
    plink_map[0] = plink_map[0].astype(str)
    plink_map[3] = plink_map[3].astype(int)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    plink_map.to_csv(args.output, sep="\t", header=False, index=False)
    print(f"PLINK .map file written to {args.output}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
