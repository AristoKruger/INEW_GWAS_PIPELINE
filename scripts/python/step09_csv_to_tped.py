import argparse
from pathlib import Path

import pandas as pd

ALLELE_MAP = {
    "0": "A A",
    "1": "A B",
    "2": "B B",
    "9": "0 0",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert a transposed genotype CSV into a PLINK TPED file",
    )
    parser.add_argument(
        "--genotypes",
        default="outputs/ibs/filtered_genotypes_strict_T.csv",
        type=Path,
        help="Transposed genotype CSV (Taxa column + SNP columns)",
    )
    parser.add_argument(
        "--map",
        default="outputs/plink/filtered_genotypes_strict.map",
        type=Path,
        help="PLINK .map file used to order SNPs",
    )
    parser.add_argument(
        "--output",
        default="outputs/plink/filtered_genotypes_strict.tped",
        type=Path,
        help="Destination TPED file",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    geno_df = pd.read_csv(args.genotypes, index_col=0)
    map_df = pd.read_csv(args.map, sep="\t", header=None)
    map_df.columns = ["CHR", "SNP", "CM", "BP"]

    missing_snps = map_df["SNP"].difference(geno_df.columns)
    if not missing_snps.empty:
        print(f"Warning: {len(missing_snps)} SNPs from map not found in genotype matrix")

    snp_order = [s for s in map_df["SNP"] if s in geno_df.columns]
    geno_df = geno_df[snp_order]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in map_df.itertuples(index=False):
            snp = row.SNP
            if snp not in geno_df.columns:
                continue
            genotypes = geno_df[snp].astype(str).map(ALLELE_MAP).fillna("0 0")
            line = f"{row.CHR}\t{snp}\t{row.CM}\t{row.BP}\t" + "\t".join(genotypes)
            handle.write(line + "\n")
    print(f"TPED written to {args.output}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
