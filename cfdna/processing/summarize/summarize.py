from ..cnv.segmentation import call_cnv_pipline
from ..fragmentation.frag_pattern import fragment_profile
from ..nucleosome.nfr import summarize_nfr
import os


def summarize(cfdna_object, frags, cnv_binsize=100000, hmm_binsize=1000000, method="bcp_online_both", outlier_smooth=True,
                  gauss_smooth=False, bcp_cutoff=0.3, nfr_filenames=["hg19_TSS.bed.gz","hg19_CTCF.bed.gz","hg19_TFBS.bed.gz"],
                  stranded=[True,False,False], normal=[0.1, 0.5, 0.9], ploidy=[2], estimatePloidy=False, minSegmentBins=25,
                  maxCN=7, verbose=False):
    """
    Calculate summary metrics

    Parameters
    ----------
        bin_size
            int
        genome
            str
        bin_bias_h5_fn
            str
        outlier_smooth
            bool
        gauss_smooth
            bool
        n_threads
            int
        tss_fn
            str
        normal
            list
        ploidy
            list
        estimatePloidy
            bool
        minSegmentBins
            int
        maxCN
            int
        bcp_cutoff
            float

    Returns
    -------
        metrics
            cfDNA_summary
    """

    # Determine file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]

    # Check if file was previously annotated
    if file_name not in cfdna_object.filenames:
        cfdna_object.log_filename(frags)

    # Call CNVs
    call_cnv_pipline(cfdna_object, frags, cnv_binsize=cnv_binsize, hmm_binsize=hmm_binsize, method=method, outlier_smooth=outlier_smooth,
                         gauss_smooth=gauss_smooth, bcp_cutoff=bcp_cutoff, normal=normal, ploidy=ploidy, estimatePloidy=estimatePloidy,
                         minSegmentBins=minSegmentBins, maxCN=maxCN, verbose=verbose)

    # Determine fragmentation profile
    fragment_profile(cfdna_object, frags, bin_size=cnv_binsize, bin_bias_h5_fn=None,
                     bias_correct=True, smooth=True, verbose=verbose)
    
    # Calculate NFR
    summarize_nfr(cfdna_object, frags, "long", nfr_filenames=nfr_filenames, stranded=stranded,
                  protection=120, min_length=120, max_length=220,
                  cnv_binsize=cnv_binsize, method="mean", peak_min_length=50,
                  peak_max_length=250, progress_bar=True, verbose=verbose)

    summarize_nfr(cfdna_object, frags, "short", nfr_filenames=nfr_filenames, stranded=stranded,
                  protection=16, min_length=16, max_length=120,
                  cnv_binsize=cnv_binsize, method="mean", peak_min_length=50,
                  peak_max_length=200, progress_bar=True, verbose=verbose)



    