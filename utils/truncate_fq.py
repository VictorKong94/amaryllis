# Run with command:
#   python3 truncate_fq.py \
#           <original_file_to_truncate> \
#           <new_truncated_file_to_create> \
#           <number_of_lines_to_keep>

import gzip
import os.path
import sys

def truncate(infile, outfile, num_lines):
    with gzip.open(infile, 'rb') as in_f, gzip.open(outfile, 'wb') as out_f:
        for i, line in enumerate(in_f):
            if i >= num_lines:
                break
            out_f.write(line)

if __name__ == '__main__':
    infile = sys.argv[1]
    outfile = sys.argv[2]
    assert not os.path.isfile(outfile), outfile + ' already exists'
    num_lines = int(sys.argv[3])
    truncate(infile, outfile, num_lines)
    