import ngsfragments as ngs
import numpy as np
import pandas as pd
from ailist import LabeledIntervalArray
import scipy
from scipy.signal import lfilter
from scipy.stats import trim_mean
from scipy.signal import periodogram, welch
from scipy.interpolate import Rbf

# Local import
from ...core import cfDNA


def wps(frags, chrom=None, protection=120, min_length=120, max_length=210, normalize=True, norm_method="mean"):
    """
    Calculate Window Protection Score (WPS)

    Parameters
    ----------
        frags : ngsfragments.fragments
            Fragments object
        chroms : str or list
            Name of chromosome
        protection : int
            Protection window
        min_length : int
            Minimum DNA fragment length
        max_length : int
            Maximum fragment length
        normalize : bool
            Whether to normalize
        norm_method : str
            Normalization method

    Returns
    -------
        wps_score : pandas.Series
            Window protection scores
    """
    
    # Caclulate WPS
    wps_score = frags.wps(chrom=chrom, protection=protection, min_length=min_length, max_length=max_length)
    if normalize:
        ngs.utilities.normalize_wps(wps_score, method=norm_method)

    return wps_score


def coverage(frags, chrom=None, min_length=120, max_length=210, normalize=True, norm_method="mean"):
    """
    Calculate Window Protection Score (WPS)

    Parameters
    ----------
        frags : ngsfragments.fragments
            Fragments object
        chroms : str or list
            Name of chromosome
        min_length : int
            Minimum DNA fragment length
        max_length : int
            Maximum fragment length
        normalize : bool
            Whether to normalize
        norm_method : str
            Normalization method

    Returns
    -------
        cov : pandas.Series
            Coverage
    """
    
    # Caclulate coverage
    cov = frags.coverage(chrom=chrom, min_length=min_length, max_length=max_length)
    if normalize:
        ngs.utilities.normalize_wps(cov, method=norm_method)

    return cov


def predict_nucleosomes(fragments,
                        protection = 120,
                        merge_distance = 5,
				        min_length = 120,
                        max_length = 500,
                        verbose = False):
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

    # Iterate over chromosomes
    chroms = fragments.frags.unique_labels

    # Initialize peaks
    total_peaks = LabeledIntervalArray()

    # Iterate over chromosomes
    for chrom in chroms:
        if verbose: print(chrom, flush=True)
        wps = fragments.frags.wps(protection, chrom, min_length, max_length)
        wps[chrom].values[:] = ngs.peak_calling.CallPeaks.normalize_signal(wps[chrom].values)

        peaks = ngs.peak_calling.CallPeaks.call_peaks(wps[chrom], str(chrom), merge_distance, 50, 190)
        if peaks.size < 10:
            continue
        midpoints = (peaks.starts + ((peaks.ends - peaks.starts) / 2)).astype(int)
        standard_peaks = LabeledIntervalArray()
        standard_peaks.add(midpoints - 84, midpoints + 84, np.repeat(chrom, len(midpoints)))
        standard_peaks = standard_peaks.merge()
        midpoints = (standard_peaks.starts + ((standard_peaks.ends - standard_peaks.starts) / 2)).astype(int)
        standard_peaks = LabeledIntervalArray()
        standard_peaks.add(midpoints - 84, midpoints + 84, np.repeat(chrom, len(midpoints)))
        total_peaks.append(standard_peaks)

    return total_peaks


def merge_nucleosomes(peaks1, peaks2, merge_distance=10):
    """
    Merge two sets of nucleosomes

    Parameters
    ----------
        peaks1
            dict[chrom:AIList]
        peaks2
            dict[chrom:AIList]
        merge_distance
            int

    Returns
    -------
        merged_peaks
            dict[chrom:AIList]
    """

    # Merge peaks
    merged_peaks = peaks1.union(peaks2).merge(merge_distance)
    
    # Re-center peaks
    midpoints = (merged_peaks.starts + ((merged_peaks.ends - merged_peaks.starts) / 2)).astype(int)
    standard_peaks = LabeledIntervalArray()
    standard_peaks.add(midpoints - 84, midpoints + 84, merged_peaks.labels)

    return standard_peaks


def internucleosomal_distances(peaks, max_distance=1000):
    """
    Determine distance between nucleosomes
    """

    # Iterate over chromosomes
    chroms = peaks.unique_labels

    # Initialize distances
    distances = []

    # Iterate over chromosomes
    for chrom in chroms:
        chrom_peaks = peaks.get(chrom)
        midpoints = (chrom_peaks.starts + ((chrom_peaks.ends - chrom_peaks.starts) / 2)).astype(int)

        pdist = np.diff(midpoints)
        distances.extend(pdist[pdist < max_distance])

    distances = np.array(distances)

    return distances


def peak_distances(self, fragments, peaks=None, bin_size=100000, max_distance=10000):
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
        peaks = self.predict_nucleosomes(fragments)

    # Create p_dist IntervalFrame
    starts = np.array([],dtype=int)
    chroms = np.array([], dtype="U25")
    for chrom in peaks:
        new_starts = np.arange(0, self.chrom_sizes[chrom]+bin_size, bin_size)
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


def spec_pgram(x, smoothed=False):
    """
    """
    
    # Compute the periodogram
    frequencies, power_spectrum = periodogram(x,
                                              window="hann",
                                              nfft=None,
                                              return_onesided=True,
                                              scaling="density",
                                              detrend="linear")
    if smoothed:
        # Smooth the periodogram using Welch's method
        frequencies, power_spectrum = welch(x,
                                          window="hann",
                                          nperseg=None,
                                          noverlap=None,
                                          detrend="linear",
                                          scaling="density",
                                          return_onesided=True)
    return frequencies, power_spectrum


def wps_gene_fft(fragments,
                 genome_version = "hg38",
                 protection = 120,
				 min_length = 120,
                 max_length = 1000,
                 freq_range = (120, 280),
                 scale = True,
                 feature = "gene",
                 verbose = False):
    """
    """

    import genome_info

    genome = genome_info.GenomeInfo(genome_version)
    #genes = genome.get_intervals("tss", downstream=10000, filter_column="gene_type", filter_selection="protein_coding")
    if feature == "gene":
        genes = genome.get_intervals("gene", filter_column="gene_type", filter_selection="protein_coding")
    elif feature == "tss":
        genes = genome.get_intervals("tss", upstream=5000, downstream=5000, filter_column="gene_type", filter_selection="protein_coding")
    elif feature == "gene_body":
        genes = genome.get_intervals("gene", upstream=5000, filter_column="gene_type", filter_selection="protein_coding")

    # Iterate over chromosomes
    chroms = fragments.frags.unique_labels

    # Initialize scores
    scores = pd.DataFrame(np.zeros((genes.shape[0], freq_range[1]-freq_range[0])), index=genes.loc[:,"gene_name"].values)
    #scores = {}
    b = np.array([1])
    a = np.append(1, -(1 / np.arange(5, 100, 4)))

    # Iterate over chromosomes
    k = 0
    for chrom in chroms:
        if verbose: print(chrom, flush=True)
        wps = fragments.frags.wps(protection, chrom, min_length, max_length)
        wps[chrom].values[:] = ngs.peak_calling.CallPeaks.normalize_signal(wps[chrom].values)

        chrom_genes = genes.loc[chrom,:]

        for i in range(chrom_genes.shape[0]):
            interval = chrom_genes.index[i]
            tdataA1 = wps[chrom].loc[interval.start:interval.end].values
            x = np.append(tdataA1[:300], tdataA1)
            tdataA1 = lfilter(b, a, x)[300:]
            tdataA1 = tdataA1 - trim_mean(tdataA1, 0.1)
            resA1 = spec_pgram(tdataA1, smoothed=False)
            freq = 1 / resA1[0][1:]
            chosen = np.logical_and(freq >= freq_range[0], freq <= freq_range[1])
            if chosen.sum() < 5:
                continue
            freq = np.round(freq[chosen])
            fft_values = pd.Series(resA1[1][1:][chosen], index=freq)
            fft_values = fft_values.groupby(by=fft_values.index.values.astype(int)).mean()
            #print(fft_values, fft_values.shape, chosen, chosen.sum(), flush=True)
            interp_func = Rbf(fft_values.index, fft_values.values, smooth=0.1)
            fft_values = pd.Series(interp_func(np.arange(freq_range[0]-5, freq_range[1]+5)), index = np.arange(freq_range[0]-5, freq_range[1]+5))
            fft_values = fft_values.rolling(window=5, center=True).mean()
            fft_values = fft_values.loc[np.arange(freq_range[0], freq_range[1])]
            if scale:
                fft_values = (fft_values - fft_values.min()) / (fft_values.max() - fft_values.min())
            scores.loc[chrom_genes.df.loc[:,"gene_name"].values[i],:] = fft_values.values

    # Edit range
    scores.columns = np.arange(freq_range[0], freq_range[1])

    return scores


def wps_windows(cfdna_object: cfDNA,
                frags: ngs.Fragments,
                key: str = "tss_wps",
                protection: int = 120,
                min_length: int = 120,
                max_length: int = 210,
                feature: str = "tss",
                smooth: bool = False) -> cfDNA:
    """
    Calculate Window Protection Score (WPS) for windows

    Parameters
    ----------
        cfdna_object : cfDNA
            cfDNA object
        frags : ngs.Fragments
            Fragments object
        key : str
            Key for storing the WPS (default: "tss_wps")
        protection : int
            Protection window (default: 120)
        min_length : int
            Minimum DNA fragment length (default: 120)
        max_length : int
            Maximum fragment length (default: 210)
        feature : str
            Feature (default: "tss")
        smooth : bool
            Whether to smooth (default: False)
        
    Returns
    -------
        cfdna_object : cfDNA
            cfDNA object

    
    """

    scores = ngs.metrics.wps_windows(frags,
                                    protection = protection,
                                    min_length = min_length,
                                    max_length = max_length,
                                    feature = feature,
                                    smooth = smooth)
    
    cfdna_object.add_obs_values(list(cfdna_object.obs)[0], key, scores)

    return cfdna_object