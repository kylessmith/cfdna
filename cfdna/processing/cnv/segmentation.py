import ngsfragments as ngs
import os
import numpy as np
from typing import List
np.seterr(all="ignore")
#from ngsfragments.segment.cnv import call_cnvs, process_cnvs
#from ngsfragments.segment import cnv_utilities

# Local imports
from ...core.core import cfDNA


def call_cnvs(cfdna_object: cfDNA,
              frags: ngs.Fragments,
              bin_size: int = 100000,
              method: str = "bcp_online_both",
              outlier_smooth: bool =True,
              gauss_smooth: bool = False,
              bcp_cutoff: float = 0.3,
              hazard: int = 100,
              shuffles: int = 5000,
              p: float = 0.00005,
              verbose: bool = False,
              **kwargs):
    """
    Call Copy Number Variants (CNVs)

    Parameters
    ----------
    chrom : str
        Name of chromosome
    bin_size : int
    min_length : int
    max_length : int
    genome : str
    genome_fn : str
    centro_fn : str
    detail_fn : str
    outlier_smooth : bool
    gauss_smooth : bool
    normal : list
    ploidy : list
    estimatePloidy : bool
    minSegmentBins : int
    maxCN : int
    bcp_cutoff : float
        
    Returns
    -------
        cnv : :class:`pandas.DataFrame`

    """

    # Find file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]

    # Log fragments
    cfdna_object.log_fragments(frags)

    # Log bin size
    cfdna_object.params["cnv_binsize"] = bin_size

    # Call Copy Number Variations
    cnv_bins, cnv_segs = ngs.segment.cnv.call_cnvs(frags,
                                                   bin_size = bin_size,
                                                   genome = str(cfdna_object.params["ref_genome"]),
                                                   method = method,
                                                   outlier_smooth = outlier_smooth,
                                                   gauss_smooth = gauss_smooth,
                                                   verbose = verbose,
                                                   cutoff = bcp_cutoff,
                                                   hazard = hazard,
                                                   shuffles = shuffles,
                                                   p = p)
    
    # Append bins
    cnv_bins = cnv_bins.loc[:,["ratios"]]
    cnv_bins.df.columns = [file_name]
    cfdna_object.add_obs_intervals(file_name, "cnv_segments", cnv_segs)
    cfdna_object.add_intervals("cnv_bins", cnv_bins)


def predict_purity(cfdna_object: cfDNA,
                   frags: ngs.Fragments,
                   file_name: str = None,
                   bin_size: int = 1000000,
                   method: str = "bcp_online_both",
                   outlier_smooth: bool = True,
                   gauss_smooth: bool = False,
                   bcp_cutoff: float = 0.3,
                   normal: List[float] = [0.1, 0.5, 0.9],
                   ploidy: List[int] = [1, 2, 3],
                   estimatePloidy: bool = False,
                   minSegmentBins: int = 25,
                   maxCN: int = 7,
                   verbose: bool = False,
                   **kwargs):
    """
    Predict tumor purity
    """

    # Find file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]

    # Check if file was previously annotated
    cfdna_object.log_fragments(frags)

    # Calculate bins
    hmm_binsize = bin_size
    bins, segments = ngs.segment.cnv.call_cnvs(frags,
                                               bin_size = bin_size,
                                               genome = str(cfdna_object.params["ref_genome"]),
                                               method = method,
                                               outlier_smooth = outlier_smooth,
                                               gauss_smooth = gauss_smooth,
                                               verbose = verbose,
                                               cutoff = bcp_cutoff)

    # Define gender
    try:
        gender = cfdna_object.obs_values[file_name]["gender"]
    except KeyError:
        gender = None
    
    # Classify calls be HMM
    hmm_states = ngs.segment.cnv_utilities.train_hmm(bins,
                                                     normal = normal,
                                                     ploidy = ploidy,
                                                     gender = gender,
                                                     estimatePloidy = estimatePloidy,
                                                     minSegmentBins = minSegmentBins,
                                                     maxCN = maxCN,
                                                     verbose = verbose,
                                                     **kwargs)

    # Assign attributes
    cfdna_object.add_anno("purity", file_name, 1 - hmm_states["n"])
    cfdna_object.add_anno("ploidy", file_name, hmm_states["phi"])
    cfdna_object.add_anno("clonal", file_name, hmm_states["Frac_genome_subclonal"])


def train_hmm(cfdna_object: cfDNA,
              frags: ngs.Fragments,
              bin_size: int = 1000000,
              method: str = "bcp_online_both",
              outlier_smooth: bool = True,
              gauss_smooth: bool = False,
              bcp_cutoff: float = 0.3,
              normal: List[float] = [0.1, 0.5, 0.9],
              ploidy: List[int] = [2],
              estimatePloidy: bool = False,
              minSegmentBins: int = 25,
              maxCN: int = 7,
              verbose: bool = False,
              **kwargs):
    """
    Train Hidden Markov Model (HMM)
    """

    # Determine file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]

    # Check if file was previously annotated
    cfdna_object.log_fragments(frags)

    # Determine bin size
    cfdna_object.params["hmm_binsize"] = bin_size

    # Calculate bins
    bins, segments = ngs.segment.cnv.call_cnvs(frags,
                                               bin_size = bin_size,
                                               genome = str(cfdna_object.params["ref_genome"]),
                                               method = method,
                                               outlier_smooth = outlier_smooth,
                                               gauss_smooth = gauss_smooth,
                                               verbose = verbose,
                                               cutoff = bcp_cutoff)

    # Define gender
    try:
        gender = cfdna_object.obs_values[file_name]["gender"]
    except KeyError:
        gender = None

    # Train HMM
    hmm_states = ngs.segment.cnv_utilities.train_hmm(bins,
                                                     normal = normal,
                                                     ploidy = ploidy,
                                                     gender = gender,
                                                     estimatePloidy = estimatePloidy,
                                                     minSegmentBins = minSegmentBins,
                                                     maxCN = maxCN,
                                                     verbose = verbose,
                                                     **kwargs)

    cfdna_object.uns[file_name] = {"hmm_states": hmm_states}



def process_cnvs(cfdna_object: cfDNA,
                 frags: ngs.Fragments = None,
                 file_name: str = None,
                 merge: bool = True):
    """
    """

    # Find file_name
    if frags is None and file_name is None:
        raise NameError("frags or file_name must be provided.")
    elif frags is not None:
        path = os.path.normpath(frags.sam_file)
        file_name = path.split(os.sep)[-1]

        # Check if file was previously annotated
        cfdna_object.log_fragments(frags)
    
    # Check if hmm has been run
    try:
        cfdna_object.uns[file_name]["hmm_states"]
    except KeyError:
        raise AttributeError("Must run .train_hmm() before .process_cnvs()")

    # Extract sample bins
    cfdna_sample = cfdna_object.intervals["cnv_bins"].loc[:,[file_name]]
    cfdna_sample.df.columns = ["ratios"]

    # Merge segments
    if merge:
        merged = ngs.segment.cnv_utilities.merge_segments(cfdna_object.obs_intervals[file_name]["cnv_segments"],
                                                          cfdna_sample)
        cfdna_object.obs_intervals[file_name]["cnv_segments"] = merged

    # Process CNVs
    processed_cnvs = ngs.segment.cnv_utilities.hmm_classify(cfdna_object.obs_intervals[file_name]["cnv_segments"],
                                                            cfdna_object.uns[file_name]["hmm_states"])
    ngs.segment.cnv_utilities.validate_distributions(processed_cnvs,
                                                     cfdna_sample)
    cfdna_object.add_obs_intervals(file_name, "cnv_segments", processed_cnvs)


def call_cnv_pipeline(cfdna_object: cfDNA,
                     frags: ngs.Fragments,
                     cnv_binsize: int = 100000,
                     hmm_binsize: int = 1000000,
                     method: str = "bcp_online_both",
                     outlier_smooth: str = True,
                     gauss_smooth: bool = False,
                     bcp_cutoff: float = 0.3,
                     normal: List[float] = [0.1, 0.25, 0.5, 0.75, 0.9],
                     ploidy: List[int] = [2],
                     estimatePloidy: bool = False,
                     minSegmentBins: int = 25,
                     maxCN: int = 5,
                     hazard: int = 100,
                     shuffles: int = 5000,
                     p: float = 0.00005,
                     verbose: bool = False):
    """
    """

    # Determine file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]

    # Check if file was previously annotated
    cfdna_object.log_fragments(frags)

    # Predict purity and ploidy
    predict_purity(cfdna_object,
                   frags,
                   bin_size = hmm_binsize,
                   method = method,
                   outlier_smooth = outlier_smooth,
                   gauss_smooth = gauss_smooth,
                   bcp_cutoff = bcp_cutoff,
                   normal = normal,
                   ploidy = ploidy,
                   estimatePloidy = estimatePloidy,
                   minSegmentBins = minSegmentBins,
                   maxCN = maxCN,
                   verbose = verbose)

    # Train HMM
    train_hmm(cfdna_object,
              frags,
              bin_size = hmm_binsize,
              method = method,
              outlier_smooth = outlier_smooth,
              gauss_smooth = gauss_smooth,
              bcp_cutoff = bcp_cutoff,
              normal = normal,
              ploidy = ploidy,
              estimatePloidy = estimatePloidy,
              minSegmentBins = minSegmentBins,
              maxCN = maxCN,
              verbose = verbose)

    # Call Copy Number Variations
    call_cnvs(cfdna_object,
              frags,
              bin_size = cnv_binsize,
              method = method,
              outlier_smooth = outlier_smooth,
              gauss_smooth = gauss_smooth,
              verbose = verbose,
              cutoff = bcp_cutoff,
              hazard = hazard,
              shuffles = shuffles,
              p = p)
    process_cnvs(cfdna_object,
                 frags)

    # Calculate segment variance
    cfdna_object.obs_intervals[file_name]["cnv_segments"].annotate(cfdna_object.intervals["cnv_bins"], file_name, "var")