from ngsfragments import fragments


def readBam(bam_filename, min_size=10, max_size=1000,
            paired=True, qcfail=False, mapq_cutoff=25,
            proportion=1.0, verbose=False):
    """
    Read BAM file

    Parameters
    ----------
        bam_filename : str
            Name of the BAM file to read
        min_size : int
            Minimum length of fragment to consider
        max_size : int
            Maximum length of fragment to consider
        paired : bool
            Whether to process as paired reads
        qcfail : bool
            Whether to keep QC failed reads
        mapq_cutoff : int
            Minimum mapping quality allowed
        proportion : float
            Proportion of reads to keep
        verbose : bool
            Print process
    
    Returns
    -------
        fragments : :class: `ngsfragments.ngsfragments`
            DNA fragment intervals from BAM file

    """

    # Read fragments from bam file
    frags = fragments(bam_filename, min_size=min_size, max_size=max_size,
                        paired=paired, qcfail=qcfail, mapq_cutoff=mapq_cutoff, verbose=verbose,
                        proportion=proportion)

    return frags