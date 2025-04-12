# First argument should be name of directory that contains the TBS files and nothing else.
# FASTQ files should be downloaded before this.

import os
import sys
assert(len(sys.argv) == 2)
dir_path = sys.argv[1]
for old_fname in os.listdir(dir_path):
    if old_fname.endswith(".fastq.gz"):
        parts = old_fname.split('_')
        new_fname = f"{parts[0]}_{parts[3]}.fastq.gz"
        #print(os.path.join(dir_path, old_fname), os.path.join(dir_path, new_fname))
        os.rename(os.path.join(dir_path, old_fname), os.path.join(dir_path, new_fname))
