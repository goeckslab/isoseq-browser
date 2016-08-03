#
# Simple script to remove reads from a SAM file stream that are not aligned to
# the provided set of chromosomes/contigs.
#

import sys

# Get chromosome names for filtering.
names = set( line.strip() for line in open( sys.argv[1] ).readlines() )

# Filter SAM file so that only reads from valid contigs are included.
for line in sys.stdin:
    if line.startswith("@"):
        print line,
    chr = line.split("\t")[2]
    if chr in names:
        print line,
