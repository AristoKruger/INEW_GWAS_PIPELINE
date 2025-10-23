from __future__ import division
import sys
import csv

if len(sys.argv) != 3:
    print("Usage: python step05_calculate_ibs.py <input_transposed_csv> <output_ibs_matrix_csv>")
    sys.exit(1)

geno_file = open(sys.argv[1], 'r')
ibs_file = open(sys.argv[2], 'w', newline='')  # Python 3 fix

# Step 1: Parse input
genotypes = []
genotype_data = {}
genotype_data_eval = {}

header = True
snp_total = 0
snp_counted = False

for line in geno_file:
    if header:
        header = False
    else:
        line2 = line.rstrip('\n').split(',')
        genotype_id = line2[0]
        genotype_snps = line2[1:]
        genotypes.append(genotype_id)
        genotype_data[genotype_id] = genotype_snps

        eval_flags = [snp in {'0', '1', '2'} for snp in genotype_snps]
        genotype_data_eval[genotype_id] = eval_flags

        if not snp_counted:
            snp_total = len(genotype_snps)
            snp_counted = True

# Step 2: Prepare results matrix
psm_results = [None] * (len(genotypes) + 1)
psm_results[0] = ['genotype'] + genotypes

# Step 3: Compute IBS matrix
for y in range(len(genotypes)):
    genotype1 = genotypes[y]
    geno1 = genotype_data[genotype1]
    geno1_eval = genotype_data_eval[genotype1]

    psm_results[y + 1] = [None] * (len(genotypes) + 1)
    psm_results[y + 1][0] = genotype1

    for x in range(y + 1):
        genotype2 = genotypes[x]
        geno2 = genotype_data[genotype2]
        geno2_eval = genotype_data_eval[genotype2]

        common = 0
        scored = 0

        for snp_number in range(snp_total):
            g1, g2 = geno1[snp_number], geno2[snp_number]
            valid1, valid2 = geno1_eval[snp_number], geno2_eval[snp_number]

            if g1 == g2 and valid1:
                common += 1
                scored += 1
            elif g1 == '1' and g2 in {'0', '2'} and valid1:
                common += 0.5
                scored += 1
            elif g2 == '1' and g1 in {'0', '2'} and valid1:
                common += 0.5
                scored += 1
            elif valid1 and valid2:
                scored += 1

        similarity = common / scored if scored else 0
        psm_results[y + 1][x + 1] = similarity
        psm_results[x + 1][y + 1] = similarity

# Step 4: Write output
writer = csv.writer(ibs_file)
writer.writerows(psm_results)
geno_file.close()
ibs_file.close()

print(f"✅ IBS matrix written to {sys.argv[2]}")

