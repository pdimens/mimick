# Features
- Swap to BGZFLib.jl and enable b/gzip decompression of input FASTA files
- rm fasta indexing
- print info text to stderr when contig length < avg mol size
- add `mask-n` to allow Ns in the input fasta to be overwritten with random ATCG bases
- intelligent insert-finder that uses reliable math to position random fragments instead of random positions with attempts

# Fixes
- None

# Breaking Changes
- None
