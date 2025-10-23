# step15_fix_pruned_ped.py
#
# Normalise the PLINK PED emitted by step11 so TASSEL receives the expected coding.
# - Increment the sex column by 1 (0 -> 1, 1 -> 2, etc.).
# - Replace allele code "B" with "C" across the genotype columns.
# Usage: python scripts/python/step15_fix_pruned_ped.py <input.ped> <output.ped>

from __future__ import annotations

import sys
from pathlib import Path


def _transform_line(parts: list[str]) -> list[str]:
    if len(parts) < 6:
        raise ValueError("PED line must have at least 6 columns (FID IID PID MID SEX PHENOTYPE).")

    try:
        sex_value = int(parts[4])
    except ValueError as exc:  # pragma: no cover - defensive guard
        raise ValueError(f"Could not parse sex column value '{parts[4]}' as integer.") from exc
    parts[4] = str(sex_value + 1)

    for idx in range(6, len(parts)):
        allele = parts[idx]
        if allele.upper() == "B":
            parts[idx] = "C"
        elif allele.upper() == "C":
            # Keep existing C as uppercase
            parts[idx] = "C"
        else:
            # Preserve original value (A, 0, -9, etc.)
            pass
    return parts


def fix_ped(input_path: Path, output_path: Path) -> None:
    if not input_path.exists():
        raise FileNotFoundError(f"Input PED file not found: {input_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", encoding="utf-8") as src, output_path.open("w", encoding="utf-8") as dst:
        for line_num, line in enumerate(src, 1):
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split()
            try:
                updated = _transform_line(parts)
            except ValueError as exc:
                raise ValueError(f"Failed processing line {line_num}: {exc}") from exc
            dst.write(" ".join(updated) + "\n")


def main(args: list[str]) -> int:
    if len(args) != 2:
        print("Usage: python scripts/python/step15_fix_pruned_ped.py <input.ped> <output.ped>")
        return 1
    src = Path(args[0])
    dst = Path(args[1])
    fix_ped(src, dst)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
