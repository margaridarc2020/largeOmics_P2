import sys
import os

input_files = sys.argv[1:-1]

output_file = sys.argv[-1]

with open(output_file, 'w') as f:
    f.write('\n'.join(input_files))

