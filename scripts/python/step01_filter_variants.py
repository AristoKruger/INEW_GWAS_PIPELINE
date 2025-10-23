import sys
import csv

if len(sys.argv) != 5:
    print("Usage: python step01_filter_variants.py <input_csv> <output_csv> <maf_cutoff_percent> <missing_data_cutoff_percent>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]
maf_cutoff = float(sys.argv[3]) / 100
md_cutoff = float(sys.argv[4]) / 100

with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.reader(infile)
    writer = csv.writer(outfile)

    header = next(reader)
    writer.writerow(header)  # write header

    for row in reader:
        snp_name = row[0]
        genotype_values = row[1:]

        # Filter out missing calls (coded as '9' or blank)
        called_genotypes = [int(g) for g in genotype_values if g in {'0', '1', '2'}]

        total_alleles = 2 * len(genotype_values)
        missing_count = sum(1 for g in genotype_values if g not in {'0', '1', '2'}) * 2
        called_count = total_alleles - missing_count

        if called_count == 0:
            print(f"[SKIPPED] {snp_name}: No calls found")
            continue

        perc_missing = missing_count / total_alleles
        if perc_missing > md_cutoff:
            print(f"[SKIPPED] {snp_name}: Missing data {perc_missing*100:.2f}% exceeds cutoff")
            continue

        # Calculate allele count
        ac = 0  # alternate allele count
        for g in genotype_values:
            if g == '0':
                ac += 2
            elif g == '1':
                ac += 1
            # '2' means homozygous alternate → 0 contribution to ac

        af = ac / called_count
        maf = af if af < 0.5 else 1 - af

        if maf >= maf_cutoff:
            writer.writerow(row)
        else:
            print(f"[SKIPPED] {snp_name}: MAF {maf:.4f} below cutoff")

print(f"\n✅ Filtering complete. Output written to: {output_file}")

