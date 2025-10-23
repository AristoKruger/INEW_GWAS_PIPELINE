# step14_convert_to_structure.py
#
# Converts 012-coded genotypes into diploid format for STRUCTURE
# Usage: python scripts/python/step14_convert_to_structure.py <input.txt> <output.txt>

import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, "r") as fin, open(outfile, "w") as fout:
    for line in fin:
        fields = line.strip().split(",")
        genotype_label = fields[0]
        genotype_data = fields[1:]

        diploid_fields = [genotype_label]
        for value in genotype_data:
            if value == "0":
                diploid_fields.extend(["1", "1"])
            elif value == "1":
                diploid_fields.extend(["1", "2"])
            elif value == "2":
                diploid_fields.extend(["2", "2"])
            elif value == "9":
                diploid_fields.extend(["-9", "-9"])
            else:
                diploid_fields.extend(["ERR", "ERR"])  # Flag unexpected values

        fout.write("\t".join(diploid_fields) + "\n")

