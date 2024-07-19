from ngsfragments import Fragments
import ngsfragments as ngs


def read_sam(bam_filename: str,
            min_size: int = 10,
            max_size: int = 1000,
            paired: bool = True,
            qcfail: bool = False,
            mapq_cutoff: int = 25,
            proportion: float = 1.0,
            n_frags: int = None,
            nthreads: int = 1,
            verbose: bool = False) -> Fragments:
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
        n_frags : int
            Number of fragments to keep
        nthreads : int
            Number of threads to use
        verbose : bool
            Print process
    
    Returns
    -------
        fragments : :class: `ngsfragments.ngsfragments`
            DNA fragment intervals from BAM file

    """

    # Read fragments from bam file
    frags = ngs.read_sam.from_sam(bam_filename,
                                    min_size = min_size,
                                    max_size = max_size,
                                    paired = paired,
                                    qcfail = qcfail,
                                    mapq_cutoff = mapq_cutoff,
                                    verbose = verbose,
                                    proportion = proportion,
                                    n_frags = n_frags,
                                    nthreads = nthreads)

    return frags