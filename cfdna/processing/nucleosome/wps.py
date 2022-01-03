def wps(frags, protection=120, min_length=120, max_length=210, normalize=True, norm_method="mean"):
    """
    Calculate Window Protection Score (WPS)

    Parameters
    ----------
        chrom : str
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
        wps : pandas.Series
            Window protection scores
    """
    
    # Caclulate WPS
    wps = frags.wps(protection=protection, min_length=min_length, max_length=max_length)
    if normalize:
        utilities.normalize_wps(wps, method=norm_method)

    return wps


def predict_nucleosomes(self, fragments, wps=None, protection=120, merge_distance=5,
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
        for chrom in self.chroms:
            wps = fragments.wps(chrom=chrom, protection=protection, min_length=min_length, max_length=max_length)
            utilities.normalize_wps(wps)
            tmp_peaks = utilities.wps_peaks(fragments, wps, merge_distance=merge_distance, min_length=min_length, max_length=max_length)
            peaks[chrom] = tmp_peaks[chrom]
    else:
        peaks = utilities.wps_peaks(fragments, wps, merge_distance=merge_distance, min_length=min_length, max_length=max_length)

    return peaks


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


import ngsfragments as ngs


def wps(frags, protection=120, min_length=120, max_length=210, normalize=True, norm_method="mean"):
    """
    Calculate Window Protection Score (WPS)

    Parameters
    ----------
        chrom : str
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
        wps : pandas.Series
            Window protection scores
    """
    
    # Caclulate WPS
    wps = frags.wps(protection=protection, min_length=min_length, max_length=max_length)
    if normalize:
        ngs.utilities.normalize_wps(wps, method=norm_method)

    return wps
