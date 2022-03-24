import gzip

def get_chrom_sizes(chrom_sizes_file):

    # Initialize
    chrom_sizes = {}

    fh = get_file_handle(chrom_sizes_file)
    for line in fh:
        chrom, size = line.strip("\n").split("\t")
        chrom_sizes.setdefault(chrom, tuple([0, size]))
    fh.close()

    return(chrom_sizes)