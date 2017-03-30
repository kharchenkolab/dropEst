#!/usr/bin/python

import sys

inp_name = sys.argv[1]
out_name = sys.argv[2] if len(sys.argv) > 2 else (inp_name + '.trimmed.bed')

with open(inp_name) as inp_f:
    with open(out_name, "w") as out_f:
        for line in inp_f:
            columns = line.split()
            columns[3] = columns[3][:columns[3].find('_')]
            out_f.write('\t'.join(columns) + '\n')
