import ngsfragments as ngs
import numpy as np
from intervalframe import IntervalFrame
import pandas as pd
import os
from ...utilities.h5_utilities import printProgressBar


def calculate_nfr(scores, score_name=None, bed_fn="hg19_TSS.bed", upstream=1000, downstream=1000, stranded=False):
    """
    Calculate Nucleosome Free Region (NFR)
    """

    ngs.utilities.window_scores(scores, bed_fn=bed_fn, upstream=upstream, downstream=downstream, stranded=stranded)


def window_mean_scores(frags, values, bed_fn="hg19_TSS.bed", upstream=1000, downstream=1000, stranded=True):
    """
    Calculate mean scores for each postion within a window

    Parameters
    ----------
        values
            numpy.ndarray
        bed_fn
            str
        upstream
            int
        downstream
            int
        stranded
            bool

    Returns
    -------
        mean_values
            numpy.ndarray
    """

    # Check if bed_fn is in data directory
    if os.path.exists(os.path.join(frags.data_dir, bed_fn)):
        bed_fn = os.path.join(frags.data_dir, bed_fn)
    else:
        if ~os.path.exists(bed_fn):
            raise FileExistsError("bed_fn does not exist!")
    
    # Calculate mean scores
    mean_values = frags.window_mean_scores(values, bed_fn, upstream=upstream, downstream=downstream, stranded=stranded)
    
    return mean_values


def window_nfr(values, bed_fn="hg19_TSS.bed", upstream=1000, downstream=1000, stranded=True, scale=True):
    """
    Calculate Nucleosome Free Region score per window

    Parameters
    ----------
        values
            numpy.ndarray
        bed_fn
            str
        upstream
            int
        downstream
            int
        stranded
            bool
        smooth
            bool
        scale
            bool
        order
            int

    Returns
    -------
        window_values
            pandas.DataFrame
    """

    # Check if bed_fn is in data directory
    if os.path.exists(os.path.join(fragments.data_dir, bed_fn)):
        bed_fn = os.path.join(fragments.data_dir, bed_fn)
    else:
        if ~os.path.exists(bed_fn):
            raise FileExistsError("bed_fn does not exist!")

    # Calculate scores
    scores = ngs.utilities.window_scores(values, bed_fn=bed_fn, upstream=upstream, downstream=downstream, stranded=stranded)
    
    # Calculate mean scores
    window_values = ngs.utilities.nfr_windows(scores, bed_fn, upstream=upstream, downstream=downstream, scale=scale)

    # Create pandas Series
    window_values = pd.Series(window_values, index=scores.index.values)
    
    return window_values


def predict_nucleosomes(cfdna_object, frags, wps=None, protection=120, merge_distance=5,
                        min_length=120, max_length=210):
    """
    Call peaks from WPS

    Parameters
    ----------
        chrom
            str (chromosomes to call for)

    Returns
    -------
        peaks
            dict[chrom:AIList] (WPS peaks)
    """

    # Calculate wps
    if wps is None:
        peaks = {}
        for chrom in cfdna_object.chroms:
            wps = frags.wps(chrom=chrom, protection=protection, min_length=min_length, max_length=max_length)
            ngs.utilities.normalize_wps(wps)
            tmp_peaks = ngs.utilities.wps_peaks(frags, wps, merge_distance=merge_distance, min_length=min_length, max_length=max_length)
            peaks[chrom] = tmp_peaks[chrom]
    else:
        peaks = ngs.utilities.wps_peaks(frags, wps, merge_distance=merge_distance, min_length=min_length, max_length=max_length)

    return peaks



def peak_distances(cfdna_object, frags, peaks=None, bin_size=100000, max_distance=10000):
    """
    Determine distance between peaks

    Parameters
    ----------
        peaks
            dict[chrom:AIList]
        bin_size
            int
        max_distance
            int
        smooth
            bool

    Returns
    -------
        peak_bins
            dict[pandas.Series]
    """

    # Calculate peaks
    if peaks is None:
        peaks = predict_nucleosomes(cfdna_object, frags)

    # Create p_dist IntervalFrame
    starts = np.array([],dtype=int)
    chroms = np.array([], dtype="U25")
    for chrom in peaks:
        new_starts = np.arange(0, cfdna_object.chrom_lengths[chrom] + bin_size, bin_size)
        starts = np.append(starts, new_starts)
        chroms = np.append(chroms, np.repeat(chrom, len(new_starts)))
    ends = starts + bin_size
    p_dist = IntervalFrame.from_array(starts, ends, labels=chroms)

    # Iterate over intervals
    p_dist.df.loc[:,"mean_dist"] = 0
    for i, interval in enumerate(p_dist.index):
        chrom = interval.label
        try:
            peaks[chrom]
            overlaps = peaks[chrom].intersect(interval.start, interval.end, label=chrom)
            if len(overlaps) > 3:
                p_dist.df.loc[i,"mean_dist"] = np.mean(overlaps.extract_ends()[:-1] - overlaps.extract_starts()[1:])
        except KeyError:
            pass

    # Remove those greater than max_distance
    p_dist = p_dist.iloc[p_dist.df.loc[:,"mean_dist"].values < max_distance, :]

    return p_dist


def summarize_nfr(cfdna_object, frags, key,
                  nfr_filenames=["hg19_TSS.bed.gz","hg19_CTCF.bed.gz","hg19_TFBS.bed.gz"],
                  stranded=[True,False,False], protection=120, min_length=120, max_length=220,
                  cnv_binsize=100000, method="mean", peak_min_length=50, peak_max_length=150,
                  progress_bar=False, verbose=False):
    """
    """

    # Determine file name
    path = os.path.normpath(frags.sam_file)
    file_name = path.split(os.sep)[-1]

    # Check if file was previously annotated
    if file_name not in cfdna_object.filenames:
        cfdna_object.log_filename(frags)

    # Initialize NFR scores
    n = {}
    mean_wps = {}
    for i, bed_fn in enumerate(nfr_filenames):
        title = bed_fn.split(".")[0]
        mean_wps[title] = np.zeros(2000)
        n[title] = 0

    # Determine NFR scores
    peak_dist = IntervalFrame()
    enrichment = {}
    for c, chrom in enumerate(cfdna_object.chroms):
        if progress_bar:
            printProgressBar(c, len(cfdna_object.chroms), "nfr summary " + key)
        if verbose: print(chrom)
        wps = frags.wps(chrom=chrom, protection=protection, min_length=min_length, max_length=max_length)
        ngs.utilities.normalize_wps(wps, method=method)

        # Calculate peaks
        peaks = ngs.utilities.wps_peaks(frags, wps, merge_distance=5, min_length=peak_min_length, max_length=peak_max_length)
        peak_dist = peak_dist.concat([peak_distances(cfdna_object, frags, peaks,
                                                    bin_size=cnv_binsize,
                                                    max_distance=1000)])

        for i, bed_fn in enumerate(nfr_filenames):
            title = bed_fn.split(".")[0]
            scores = ngs.utilities.window_scores(wps, bed_fn=bed_fn, upstream=1000, downstream=1000, stranded=stranded[i], verbose=verbose)
            scores = scores.groupby(by=scores.index.values).mean()
            window_values = ngs.utilities.nfr_windows(scores, scale=True)
            window_values = pd.Series(window_values, index=scores.index.values)

            try:
                enrichment[title]
            except KeyError:
                iterables = [[chrom], window_values.index.values]
                index = pd.MultiIndex.from_product(iterables)
                enrichment[title] = pd.DataFrame(window_values.values, index=index, columns=[file_name])
            try:
                enrichment[title].loc[chrom, :]
                enrichment[title].loc[chrom, file_name] = window_values
            except KeyError:
                iterables = [[chrom], window_values.index.values]
                index = pd.MultiIndex.from_product(iterables)
                new_df = pd.DataFrame(window_values.values, index=index, columns=[file_name])
                enrichment[title] = enrichment[title].append(new_df)
            
            mean_wps[title] += scores.values.sum(axis=0)
            n[title] += scores.shape[0]

    # Calculate average
    nfr_score = {}
    for title in mean_wps:
        mean_wps[title] /= n[title]
        nfr_score[title] = ngs.utilities.nfr(mean_wps[title], center_x=1000, scale=True, order=35)
    
    # Record all values
    cfdna_object.add_intervals(key + "_peak_dist", file_name, peak_dist)
    for title in mean_wps:
        cfdna_object.add_anno(key + "_" + title + "_nfr_score", file_name, nfr_score[title][0])
        cfdna_object.add_obs_values(key + "_" + title + "_mean_wps", file_name, mean_wps[title])
        cfdna_object.add_obs_values(key + "_" + title + "_nfr_enrichment", file_name, enrichment[title].T)