import ngsfragments as ngs
import os
import numpy as np
np.seterr(all="ignore")
#from ngsfragments.segment.cnv import call_cnvs, process_cnvs
#from ngsfragments.segment import cnv_utilities


def call_cnvs(cfdna_object, frags, bin_size=100000, bin_bias_h5_fn=None, method="bcp_online_both",
                  outlier_smooth=True, gauss_smooth=False, bcp_cutoff=0.3, verbose=False, **kwargs):
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
    bin_bias_h5_fn : str
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

    # Check if file was previously annotated
    if file_name not in cfdna_object.filenames:
        cfdna_object.log_filename(frags)

    # Determine bin size
    if cfdna_object.cnv_binsize is None:
        cfdna_object.cnv_binsize = bin_size
    elif cfdna_object.cnv_binsize != bin_size:
        raise AttributeError("Provided bin_size does not match previous bin_size")

    # Call Copy Number Variations
    cnv_bins, cnv_segs = ngs.segment.cnv.call_cnvs(frags, bin_size=bin_size, bin_bias_h5_fn=bin_bias_h5_fn,
                                                    genome=cfdna_object.ref_genome, method=method, outlier_smooth=outlier_smooth,
                                                    gauss_smooth=gauss_smooth, verbose=verbose, cutoff=bcp_cutoff)
    
    # Append bins
    cnv_bins = cnv_bins.loc[:,["ratios"]]
    cnv_bins.df.columns = [file_name]
    cfdna_object.add_obs_intervals("cnv_segments", file_name, cnv_segs)
    cfdna_object.add_intervals("cnv_bins", file_name, cnv_bins)

    #if cfdna_object.cnv_bins is None:
    #    cfdna_object.cnv_bins = cnv_bins
    #else:
    #    cnv_bins = cnv_bins.exact_match(cfdna_object.cnv_bins)
    #    cfdna_object.cnv_bins = cfdna_object.cnv_bins.exact_match(cnv_bins)
    #    cfdna_object.cnv_bins.df.loc[:,file_name] = cnv_bins.loc[:,file_name].values


def predict_purity(cfdna_object, frags, bin_size = 1000000, bin_bias_h5_fn=None, method="bcp_online_both",
                    outlier_smooth=True, gauss_smooth=False, bcp_cutoff=0.3, normal = [0.1, 0.5, 0.9],
                    ploidy = [1, 2, 3], estimatePloidy=False, minSegmentBins=25, maxCN=7, verbose=False, **kwargs):
    """
    Predict tumor purity
    """

    # Find file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]

    # Check if file was previously annotated
    if file_name not in cfdna_object.filenames:
        cfdna_object.log_filename(frags)

    # Calculate bins
    hmm_binsize = bin_size
    bins, segments = ngs.segment.cnv.call_cnvs(frags, bin_size=bin_size, bin_bias_h5_fn=bin_bias_h5_fn, genome=cfdna_object.ref_genome, method=method,
                                                outlier_smooth=outlier_smooth, gauss_smooth=gauss_smooth, verbose=verbose, cutoff=bcp_cutoff)

    # Define gender
    try:
        gender = cfdna_object.gender[file_name]
    except KeyError:
        gender = None
    
    # Classify calls be HMM
    hmm_states = ngs.segment.cnv_utilities.train_hmm(bins, normal=normal, ploidy=ploidy,
                                                    gender=gender, estimatePloidy=estimatePloidy,
                                                    minSegmentBins=minSegmentBins, maxCN=maxCN,
                                                    verbose=verbose, **kwargs)

    # Assign attributes
    cfdna_object.add_anno("purity", file_name, 1 - hmm_states["n"])
    cfdna_object.add_anno("ploidy", file_name, hmm_states["phi"])
    cfdna_object.add_anno("clonal", file_name, hmm_states["Frac_genome_subclonal"])
    #cfdna_object.purity[file_name] = 1 - hmm_states["n"]
    #cfdna_object.ploidy[file_name] = hmm_states["phi"]
    #cfdna_object.clonal[file_name] = hmm_states["Frac_genome_subclonal"]


def train_hmm(cfdna_object, frags, bin_size = 1000000, bin_bias_h5_fn=None, method="bcp_online_both",
                outlier_smooth=True, gauss_smooth=False, bcp_cutoff=0.3, normal = [0.1, 0.5, 0.9],
                ploidy = [2], estimatePloidy=False, minSegmentBins=25, maxCN=7, verbose=False, **kwargs):
    """
    Train Hidden Markov Model (HMM)
    """

    # Determine file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]

    # Check if file was previously annotated
    if file_name not in cfdna_object.filenames:
        cfdna_object.log_filename(frags)

    # Determine bin size
    if cfdna_object.hmm_binsize is None:
        cfdna_object.hmm_binsize = bin_size
    elif cfdna_object.hmm_binsize != bin_size:
        raise AttributeError("Provided bin_size does not match previous bin_size")

    # Calculate bins
    bins, segments = ngs.segment.cnv.call_cnvs(frags, bin_size=bin_size, bin_bias_h5_fn=bin_bias_h5_fn, genome=cfdna_object.ref_genome, method=method,
                                outlier_smooth=outlier_smooth, gauss_smooth=gauss_smooth, verbose=verbose, cutoff=bcp_cutoff)

    # Define gender
    try:
        gender = cfdna_object.gender[file_name]
    except KeyError:
        gender = None

    # Train HMM
    cfdna_object.hmm_states[file_name] = ngs.segment.cnv_utilities.train_hmm(bins, normal=normal, ploidy=ploidy,
                                                                gender=gender, estimatePloidy=estimatePloidy,
                                                                minSegmentBins=minSegmentBins, maxCN=maxCN,
                                                                verbose=verbose, **kwargs)



def process_cnvs(cfdna_object, frags=None, file_name=None, merge=True):
    """
    """

    # Find file_name
    if frags is None and file_name is None:
        raise NameError("frags or file_name must be provided.")
    elif frags is not None:
        path = os.path.normpath(frags.sam_file)
        file_name = path.split(os.sep)[-1]

        # Check if file was previously annotated
        if file_name not in cfdna_object.filenames:
            cfdna_object.log_filename(frags)
    
    # Check if hmm has been run
    try:
        cfdna_object.hmm_states[file_name]
    except KeyError:
        raise AttributeError("Must run .train_hmm() before .process_cnvs()")

    # Extract sample bins
    cfdna_sample = cfdna_object.intervals["cnv_bins"].loc[:,[file_name]]
    cfdna_sample.df.columns = ["ratios"]

    # Merge segments
    if merge:
        cfdna_object.obs_intervals["cnv_segments"][file_name] = ngs.segment.cnv_utilities.merge_segments(cfdna_object.obs_intervals["cnv_segments"][file_name],
                                                                                        cfdna_sample)

    # Process CNVs
    cfdna_object.obs_intervals["cnv_segments"][file_name] = ngs.segment.cnv_utilities.hmm_classify(cfdna_object.obs_intervals["cnv_segments"][file_name],
                                                                                  cfdna_object.hmm_states[file_name])
    ngs.segment.cnv_utilities.validate_distributions(cfdna_object.obs_intervals["cnv_segments"][file_name],
                                                     cfdna_sample)


def call_cnv_pipline(cfdna_object, frags, cnv_binsize=100000, hmm_binsize=1000000, method="bcp_online_both", outlier_smooth=True,
                         gauss_smooth=False, bcp_cutoff=0.3, normal=[0.1, 0.5, 0.9], ploidy=[2], estimatePloidy=False, minSegmentBins=25,
                         maxCN=7, verbose=False):
    """
    """

    # Determine file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]

    # Check if file was previously annotated
    if file_name not in cfdna_object.filenames:
        cfdna_object.log_filename(frags)

    # Predict purity and ploidy
    predict_purity(cfdna_object, frags, bin_size = hmm_binsize, bin_bias_h5_fn=None, method=method,
                        outlier_smooth=outlier_smooth, gauss_smooth=gauss_smooth, bcp_cutoff=bcp_cutoff, normal=normal,
                        ploidy=ploidy, estimatePloidy=estimatePloidy, minSegmentBins=minSegmentBins, maxCN=maxCN, verbose=verbose)

    # Train HMM
    train_hmm(cfdna_object, frags, bin_size = hmm_binsize, bin_bias_h5_fn=None, method=method,
                outlier_smooth=outlier_smooth, gauss_smooth=gauss_smooth, bcp_cutoff=bcp_cutoff, normal=normal,
                ploidy=ploidy, estimatePloidy=estimatePloidy, minSegmentBins=minSegmentBins, maxCN=maxCN, verbose=verbose)

    # Call Copy Number Variations
    call_cnvs(cfdna_object, frags, bin_size=cnv_binsize, bin_bias_h5_fn=None, method=method,
                outlier_smooth=outlier_smooth, gauss_smooth=gauss_smooth, verbose=verbose, cutoff=bcp_cutoff)
    process_cnvs(cfdna_object, frags)

    # Calculate segment variance
    cfdna_object.obs_intervals["cnv_segments"][file_name].annotate(cfdna_object.intervals["cnv_bins"], file_name, "var")