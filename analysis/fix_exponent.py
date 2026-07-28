#!/usr/bin/env python3

import re

input_file = "validity_all.dat"
output_file = "validity_all_fixed.dat"

# Matches numbers like:
# 0.2661808633-107
# 1.2345678900+105
# and converts them to
# 0.2661808633E-107
# 1.2345678900E+105

pat = re.compile(r'^([+-]?\d*\.?\d+)([+-]\d{2,3})$')

nfixed = 0

with open(input_file, "r") as fin, open(output_file, "w") as fout:

    for line in fin:

        # Copy header/comments/blank lines unchanged
        if line.startswith("#") or not line.strip():
            fout.write(line)
            continue

        cols = line.split()

        # Check all numerical output columns
        # Columns:
        # 0 lambda
        # 1 eta
        # 2 DG
        # 3 Exact
        # 4 HT
        # 5 ST
        # 6 Marcus
        # 7 Exact/Marcus
        # 8 HT/Marcus
        # 9 ST/Marcus
        # 10 ln(Exact/Marcus)
        # 11 ln(HT/Marcus)
        # 12 ln(ST/Marcus)

        for i in range(3, len(cols)):

            m = pat.match(cols[i])

            if m and ("E" not in cols[i]) and ("e" not in cols[i]):
                cols[i] = m.group(1) + "E" + m.group(2)
                nfixed += 1

        fout.write("    ".join(cols) + "\n")

print("Done.")
print(f"Fixed {nfixed} malformed numbers.")
print(f"Output written to: {output_file}")
