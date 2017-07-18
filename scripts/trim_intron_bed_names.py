#!/usr/bin/python

import sys

if len(sys.argv) == 1 or any([arg == '--help' or arg == '-h' for arg in sys.argv]):
    print("Usage: " + sys.argv[0] + " input_file.bed [output_file.bed]")
    exit(1)

inp_name = sys.argv[1]
out_name = sys.argv[2] if len(sys.argv) > 2 else (inp_name + '.trimmed.bed')

with open(inp_name) as inp_f:
    with open(out_name, "w") as out_f:
        for line in inp_f:
            columns = line.split()
            columns[3] = columns[3][:columns[3].find('_intron')]
            out_f.write('\t'.join(columns) + '\n')