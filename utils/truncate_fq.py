from __future__ import absolute_import, division, print_function
import gzip
import sys

def truncate(infile, num_lines):
    outfile = infile.replace('original/', 'fq/raw/')
    with gzip.open(infile, 'rb') as in_f, gzip.open(outfile, 'wb') as out_f:
        for i, line in enumerate(in_f):
            if i >= num_lines:
                break
            out_f.write(line)

if __name__ == '__main__':
    infile = sys.argv[1]
    num_lines = int(sys.argv[2])
    truncate(infile, num_lines)
