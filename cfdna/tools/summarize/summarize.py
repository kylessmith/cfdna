from ..cnv.segmentation import call_cnvs
from ..fragmentation.frag_pattern import fragment_profile
from ..nucleosome.nfr import summarize_nfr
from ...core.core import cfDNA
from ngsfragments import fragments
import os
from typing import Dict, List, Tuple


def summarize(cfdna_object: cfDNA,
              frags: fragments,
              cnv_binsize: int = 100000,
              hmm_binsize: int = 1000000,
              method: str = "bcp_online_both",
              outlier_smooth: bool = True,
              gauss_smooth: bool = False,
              bcp_cutoff: float = 0.3,
              nfr_filenames: List[str] = ["hg19_TSS.bed.gz","hg19_CTCF.bed.gz","hg19_TFBS.bed.gz"],
              stranded: List[bool] = [True,False,False],
              normal: List[float] = [0.1, 0.5, 0.9],
              ploidy: List[int] = [2],
              estimatePloidy: bool = False,
              minSegmentBins: int = 25,
              maxCN: int = 7,
              progress_bar: bool = False,
              verbose: bool = False) -> None:
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
        None: None
            None
    """

    # Determine file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]

    # Check if file was previously annotated
    if file_name not in cfdna_object.filenames:
        cfdna_object.log_filename(frags)

    # Call CNVs
    if verbose: print("Running CNV pipeline...", flush=True)
    cfdna_object = call_cnvs(cfdna_object, frags, cnv_binsize=cnv_binsize, hmm_binsize=hmm_binsize, method=method, outlier_smooth=outlier_smooth,
                            gauss_smooth=gauss_smooth, bcp_cutoff=bcp_cutoff, normal=normal, ploidy=ploidy, estimatePloidy=estimatePloidy,
                            minSegmentBins=minSegmentBins, maxCN=maxCN, verbose=verbose)

    # Determine fragmentation profile
    if verbose: print("Calculating fragment profile...", flush=True)
    fragment_profile(cfdna_object, frags, bin_size=cnv_binsize, bin_bias_h5_fn=None,
                     bias_correct=True, smooth=True, verbose=verbose)
    
    # Calculate NFR
    if verbose: print("Calling NFR for long fragments...", flush=True)
    summarize_nfr(cfdna_object, frags, "long", nfr_filenames=nfr_filenames, stranded=stranded,
                  protection=120, min_length=120, max_length=220,
                  cnv_binsize=cnv_binsize, method="mean", peak_min_length=50,
                  peak_max_length=250, progress_bar=progress_bar, verbose=verbose)

    if verbose: print("Calling NFR for short fragments...", flush=True)
    summarize_nfr(cfdna_object, frags, "short", nfr_filenames=nfr_filenames, stranded=stranded,
                  protection=16, min_length=16, max_length=120,
                  cnv_binsize=cnv_binsize, method="mean", peak_min_length=50,
                  peak_max_length=200, progress_bar=progress_bar, verbose=verbose)

    return None

    