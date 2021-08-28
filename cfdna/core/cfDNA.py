import os
import numpy as np
import pandas as pd
import h5py
import argparse
from intervalframe import IntervalFrame
from intervalframe.write.write_h5 import write_h5_intervalframe
from intervalframe.read.read_h5 import read_h5_intervalframe
from ngsfragments import fragments, utilities
from ngsfragments.segment.cnv import call_cnvs, process_cnvs
from ngsfragments.segment import cnv_utilities
from ngsfragments.segment import correction
from .. import Plot


class cfDNA(object):
    """Wrapper of fragments for cell-free DNA

	:class:`~cfDNA.cfDNA` stores a fragments object

    """

    def __init__(self, sam_fn=None, min_size=10, max_size=1000,
                 paired=True, qcfail=False, mapq_cutoff=25, proportion=1.0,
                 ref_genome="hg19", verbose=False):
        """
        Initialize cfDNA object
                 
        Params
        ------
            sam_fn : str (optional)
                 Name of the SAM file
            chrom : str
                 Name of chromosome
            sam_file : str?
                 ?
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
            n_jobs : int
                 Number of cores to use while reading sam file
                 
        Returns
        -------
            None
        """

        # Back object by with h5 file

        # Initialize analysis properties
        self.ref_genome = ref_genome
        self.gender = None
        self.l_nfr_enrichment = {}
        self.s_nfr_enrichment = {}
        self.mean_l_wps = {}
        self.mean_s_wps = {}
        self.binned_frag_profile = None
        self.binned_wps_dist = None
        self.cnv_segs = None
        self.hmm_states = None
        self.cnv_bins = None
        self.length_dist = None
        self.cnv_hmm_purity = None
        self.cnv_hmm_ploidy = None
        self.cnv_binsize = None
        self.wps_binsize = None
        self.profile_binsize = None

        # Initialize setup propteries
        self.verbose = verbose

        # Read fragments from bam file
        self.frags = fragments(sam_fn, min_size=min_size, max_size=max_size,
				          paired=paired, qcfail=qcfail, mapq_cutoff=mapq_cutoff, verbose=verbose,
				          proportion=proportion)

        # Calculate length distributions
        self.length_dist = self.frags.length_dist()


    @property
    def size(self):
        """int: Number of fragments
        """

        return self.frags.size

    @property
    def n_chroms(self):
        """int: Number of chromosomes
        """

        return self.frags.n_chroms

    @property
    def chroms(self):
        """"obj:`list` : List of chromosomes
        """

        return self.frags.chroms


    def __len__(self):
        """
        Get length from self.frags
        """

        return self.size

    
    def downsample(self, proportion=1.0, n_frags=None):
        """Randomly downsample fragments

        Parameters
        ----------
        proportion : int
            Proportion of fragments to use
        n_frags : int
            Number of fragments to use

        Returns
        -------
        :class:`~cfDNA.cfDNA`
            Contains downsamples fragments

        """

        # Initialize new cfDNA object
        filtered_cfDNA = cfDNA()

        # Filter fragments
        filtered_cfDNA.frags = self.frags.downsample(n_frags, proportion)

        # Calculate length distributions
        filtered_cfDNA.length_dist = filtered_cfDNA.frags.length_dist()

        # Reset analysis properties
        self.gender = None
        self.nfr_enrichment = None
        self.binned_frag_profile = None
        self.binned_wps_dist = None
        self.cnv_segs = None
        self.bins = None
        self.cnv_median_var = None
        self.cnv_hmm_purity = None
        self.cnv_hmm_ploidy = None
        self.cnv_binsize = None
        self.wps_binsize = None
        self.profile_binsize = None

        return filtered_cfDNA


    def wps(self, protection=120, min_length=120, max_length=210, normalize=True, norm_method="mean"):
        """
        Calculate Window Protection Score (WPS)

        Parameters
        ----------
            chrom
                str
            protection
                int
            min_length
                int
            max_length
                int
            normalize
                bool
            n_threads
                int

        Returns
        -------
            wps
                pandas.Series
        """

        wps = self.frags.wps(protection=protection, min_length=min_length, max_length=max_length)
        if normalize:
            utilities.normalize_wps(wps, method=norm_method)

        return wps


    def calculate_nfr(self, bed_fn="hg19_TSS.bed", upstream=1000, downstream=1000, stranded=False):
        """
        """

        window_scores(scores, bed_fn="hg19_TSS.bed", upstream=1000, downstream=1000, stranded=False)


    def predict_nucleosomes(self, wps=None, protection=120, merge_distance=5,
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
                wps = self.frags.wps(chrom=chrom, protection=protection, min_length=min_length, max_length=max_length)
                utilities.normalize_wps(wps)
                tmp_peaks = utilities.wps_peaks(self.frags, wps, merge_distance=merge_distance, min_length=min_length, max_length=max_length)
                peaks[chrom] = tmp_peaks[chrom]
        else:
            peaks = utilities.wps_peaks(self.frags, wps, merge_distance=merge_distance, min_length=min_length, max_length=max_length)

        return peaks

    
    def peak_distances(self, peaks=None, bin_size=100000, max_distance=10000):
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
            peaks = self.predict_nucleosomes()

        # Create p_dist IntervalFrame
        starts = np.array([],dtype=int)
        chroms = np.array([], dtype="U25")
        for chrom in peaks:
            new_starts = np.arange(0, self.frags.genome[chrom]+bin_size, bin_size)
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


    def fragment_profile(self, bin_size=100000, bin_bias_h5_fn=None, bias_correct=True, smooth=True, verbose=False):
        """
        Determine fragment profile (ratio of large fragments to small fragments)

        Parameters
        ----------
            frags : ngsfragments.fragments
            chrom
                str
            bins_size
                int
            genome
                str
            bin_bias_h5_fn
                str
            bias_correct
                bool
            smooth
                bool

        Returns
        -------
            fragment_profile
                pandas.DataFrame
        """

        # Calculate coverage of small and large fragments
        small_bin_coverage = self.frags.bin_counts(bin_size=bin_size, min_length=100, max_length=150)
        large_bin_coverage = self.frags.bin_counts(bin_size=bin_size, min_length=151, max_length=220)

        # Correct coverage
        small_bin_coverage = correction.correct(small_bin_coverage, genome=self.ref_genome, bin_size=bin_size, verbose=verbose)
        large_bin_coverage = correction.correct(large_bin_coverage, genome=self.ref_genome, bin_size=bin_size, verbose=verbose)

        # Calculate fragment profile
        fragment_profile = small_bin_coverage.copy()
        fragment_profile.df.loc[:,"counts"] = small_bin_coverage.df.loc[:,"counts"].values / large_bin_coverage.df.loc[:,"counts"].values

        # Correct regions with 0 for both
        fragment_profile.df.loc[pd.isnull(fragment_profile.loc[:,"counts"].values), "counts"] = 1.0

        return fragment_profile


    def predict_purity(self, bin_size = 1000000, bin_bias_h5_fn=None, method="bcp_online_both",
                       outlier_smooth=True, gauss_smooth=False, bcp_cutoff=0.3, normal = [0.1, 0.5, 0.9],
                       ploidy = [1, 2, 3], estimatePloidy=False, minSegmentBins=25, maxCN=7, verbose=False, **kwargs):
        """
        """

        # Calculate bins
        hmm_binsize = bin_size
        bins, segments = call_cnvs(self.frags, bin_size=bin_size, bin_bias_h5_fn=bin_bias_h5_fn, genome=self.ref_genome, method=method,
                                   outlier_smooth=outlier_smooth, gauss_smooth=gauss_smooth, verbose=verbose, cutoff=bcp_cutoff)

        # Classify calls be HMM
        hmm_states = cnv_utilities.train_hmm(bins, normal=normal, ploidy=ploidy,
                                     gender=self.gender, estimatePloidy=estimatePloidy,
                                     minSegmentBins=minSegmentBins, maxCN=maxCN,
                                     verbose=verbose, **kwargs)

        # Assign attributes
        self.purity = 1 - hmm_states["n"]
        self.ploidy = hmm_states["phi"]


    def train_hmm(self, bin_size = 1000000, bin_bias_h5_fn=None, method="bcp_online_both",
                  outlier_smooth=True, gauss_smooth=False, bcp_cutoff=0.3, normal = [0.1, 0.5, 0.9],
                  ploidy = [2], estimatePloidy=False, minSegmentBins=25, maxCN=7, verbose=False, **kwargs):
        """
        """

        # Calculate bins
        self.hmm_binsize = bin_size
        bins, segments = call_cnvs(self.frags, bin_size=bin_size, bin_bias_h5_fn=bin_bias_h5_fn, genome=self.ref_genome, method=method,
                                   outlier_smooth=outlier_smooth, gauss_smooth=gauss_smooth, verbose=verbose, cutoff=bcp_cutoff)

        # Classify calls be HMM
        self.hmm_states = cnv_utilities.train_hmm(bins, normal=normal, ploidy=ploidy,
                                     gender=self.gender, estimatePloidy=estimatePloidy,
                                     minSegmentBins=minSegmentBins, maxCN=maxCN,
                                     verbose=verbose, **kwargs)


    def call_cnvs(self, bin_size=100000, bin_bias_h5_fn=None, method="bcp_online_both",
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

        # Call Copy Number Variations
        self.cnv_binsize = bin_size
        self.cnv_bins, self.cnv_segments = call_cnvs(self.frags, bin_size=bin_size, bin_bias_h5_fn=bin_bias_h5_fn, genome=self.ref_genome, method=method,
                                                     outlier_smooth=outlier_smooth, gauss_smooth=gauss_smooth, verbose=verbose, cutoff=bcp_cutoff)
    
    
    def process_cnvs(self, merge=True):
        """
        """

        # Check if hmm has been run
        if self.hmm_states is None:
            raise AttributeError("Must run .train_hmm() before .process_cnvs()")

        # Merge segments
        if merge:
            self.cnv_segments = cnv_utilities.merge_segments(self.cnv_segments, self.cnv_bins)

        # Process CNVs
        self.cnv_segments = cnv_utilities.hmm_classify(self.cnv_segments, self.hmm_states)
        cnv_utilities.validate_distributions(self.cnv_segments, self.cnv_bins)


    def window_mean_scores(self, values, bed_fn="hg19_TSS.bed", upstream=1000, downstream=1000, stranded=True):
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
        if os.path.exists(os.path.join(self.frags.data_dir, bed_fn)):
            bed_fn = os.path.join(self.frags.data_dir, bed_fn)
        else:
            if ~os.path.exists(bed_fn):
                raise FileExistsError("bed_fn does not exist!")
        
        # Calculate mean scores
        mean_values = self.frags.window_mean_scores(values, bed_fn, upstream=upstream, downstream=downstream, stranded=stranded)
        
        return mean_values

    
    def window_nfr(self, values, bed_fn="hg19_TSS.bed", upstream=1000, downstream=1000, stranded=True, scale=True):
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
        if os.path.exists(os.path.join(self.frags.data_dir, bed_fn)):
            bed_fn = os.path.join(self.frags.data_dir, bed_fn)
        else:
            if ~os.path.exists(bed_fn):
                raise FileExistsError("bed_fn does not exist!")

        # Calculate scores
        scores = utilities.window_scores(values, bed_fn=bed_fn, upstream=upstream, downstream=downstream, stranded=stranded)
        
        # Calculate mean scores
        window_values = utilities.nfr_windows(scores, bed_fn, upstream=upstream, downstream=downstream, scale=scale)

        # Create pandas Series
        window_values = pd.Series(window_values, index=scores.index.values)
        
        return window_values


    def summarize(self, cnv_binsize=100000, hmm_binsize=1000000, method="bcp_online_both", outlier_smooth=True,
                  gauss_smooth=False, bcp_cutoff=0.3, nfr_filenames=["hg19_TSS.bed","hg19_CTCF.bed","hg19_TFBS.bed"],
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

        # Predict purity and ploidy
        self.predict_purity(bin_size = hmm_binsize, bin_bias_h5_fn=None, method=method,
                            outlier_smooth=outlier_smooth, gauss_smooth=gauss_smooth, bcp_cutoff=bcp_cutoff, normal=normal,
                            ploidy=ploidy, estimatePloidy=estimatePloidy, minSegmentBins=minSegmentBins, maxCN=maxCN, verbose=verbose)

        # Train HMM
        self.train_hmm(bin_size = hmm_binsize, bin_bias_h5_fn=None, method=method,
                  outlier_smooth=outlier_smooth, gauss_smooth=gauss_smooth, bcp_cutoff=bcp_cutoff, normal=normal,
                  ploidy=ploidy, estimatePloidy=estimatePloidy, minSegmentBins=minSegmentBins, maxCN=maxCN, verbose=verbose)

        # Call Copy Number Variations
        self.cnv_binsize = cnv_binsize
        self.cnv_bins, self.cnv_segments = call_cnvs(self.frags, bin_size=cnv_binsize, bin_bias_h5_fn=None, method=method,
                                                     outlier_smooth=outlier_smooth, gauss_smooth=gauss_smooth, verbose=verbose, cutoff=bcp_cutoff)
        self.process_cnvs()

        # Determine fragmentation profile
        self.frag_profile = self.fragment_profile(bin_size=cnv_binsize, bin_bias_h5_fn=None,
                                             bias_correct=True, smooth=True, verbose=verbose)
        
        # Initialize NFR scores
        n = {}
        for i, bed_fn in enumerate(nfr_filenames):
            title = bed_fn.split(".")[0]
            self.l_nfr_enrichment[title] = {}
            self.mean_l_wps[title] = np.zeros(2000)
            n[title] = 0
        # Determine NFR scores
        for chrom in self.chroms:
            wps = self.frags.wps(chrom=chrom, protection=120, min_length=120, max_length=220)
            utilities.normalize_wps(wps, method="mean")
            #peaks = utilities.wps_peaks(self.frags, wps, merge_distance=5, min_length=50, max_length=150)
            #self.l_peak_dist[chrom] = self.peak_distances(peaks, bin_size=cnv_binsize, max_distance=1000)
            for i, bed_fn in enumerate(nfr_filenames):
                title = bed_fn.split(".")[0]
                scores = utilities.window_scores(wps, bed_fn=bed_fn, upstream=1000, downstream=1000, stranded=stranded[i])
                scores = scores.groupby(by=scores.index.values).mean()
                window_values = utilities.nfr_windows(scores, scale=True)
                window_values = pd.Series(window_values, index=scores.index.values)
                self.l_nfr_enrichment[title][chrom] = window_values
                self.mean_l_wps[title] += scores.values.sum(axis=0)
                n[title] += scores.shape[0]

        # Calculate average
        for title in self.mean_l_wps:
            self.mean_l_wps[title] /= n[title]

        # Initialize NFR scores
        n = {}
        for i, bed_fn in enumerate(nfr_filenames):
            title = bed_fn.split(".")[0]
            self.s_nfr_enrichment[title] = {}
            self.mean_s_wps[title] = np.zeros(2000)
            n[title] = 0
        # Determine NFR scores
        for chrom in self.chroms:
            wps = self.frags.wps(chrom=chrom, protection=16, min_length=16, max_length=120)
            utilities.normalize_wps(wps, method="mean")
            for i, bed_fn in enumerate(nfr_filenames):
                title = bed_fn.split(".")[0]
                scores = utilities.window_scores(wps, bed_fn=bed_fn, upstream=1000, downstream=1000, stranded=stranded[i])
                scores = scores.groupby(by=scores.index.values).mean()
                window_values = utilities.nfr_windows(scores, scale=True)
                window_values = pd.Series(window_values, index=scores.index.values)
                self.s_nfr_enrichment[title][chrom] = window_values
                self.mean_s_wps[title] += scores.values.sum(axis=0)
                n[title] += scores.shape[0]

        # Calculate average
        for title in self.mean_s_wps:
            self.mean_s_wps[title] /= n[title]

        # 


    def write(self, h5_filename):
        """
        """
        
        # Open file
        f = h5py.File(h5_filename, "w")

        # Write attributes
        f["purity"] = self.purity
        f["ploidy"] = self.ploidy
        f["cnv_binsize"] = self.cnv_binsize
        f["wps_binsize"] = self.wps_binsize


        # Write fragment profile
        frag_profile = f.create_group("frag_profile")
        write_h5_intervalframe(self.frag_profile, frag_profile)

        # Write l nfr enrichment values
        nfr_l = f.create_group("l_nfr_enrichment")
        for title in self.l_nfr_enrichment:
            title_group = nfr_l.create_group(title)
            for chrom in self.l_nfr_enrichment[title]:
                title_group_chrom = title_group.create_group(chrom)
                title_group_chrom["nfr"] = self.l_nfr_enrichment[title][chrom].values.astype(np.float16)
                title_group_chrom["index"] = self.l_nfr_enrichment[title][chrom].index.values.astype(str).astype(bytes)

        # Write l wps mean values
        mean_l_wps = f.create_group("mean_l_wps")
        for title in self.mean_l_wps:
            mean_l_wps[title] = self.mean_l_wps[title]
        
        # Write s nfr enrichment values
        nfr_s = f.create_group("s_nfr_enrichment")
        for title in self.s_nfr_enrichment:
            title_group = nfr_s.create_group(title)
            for chrom in self.s_nfr_enrichment[title]:
                title_group_chrom = title_group.create_group(chrom)
                title_group_chrom["nfr"] = self.s_nfr_enrichment[title][chrom].values.astype(np.float16)
                title_group_chrom["index"] = self.s_nfr_enrichment[title][chrom].index.values.astype(str).astype(bytes)

         # Write s wps mean values
        mean_s_wps = f.create_group("mean_s_wps")
        for title in self.mean_s_wps:
            mean_s_wps[title] = self.mean_s_wps[title]

        # Close file
        f.close()


    def read(self, h5_filename):
        """
        """

        # Open file
        f = h5py.File(h5_filename, "r")

        # Read attributes
        purity = float(np.array(f["purity"]))
        ploidy = float(np.array(f["ploidy"]))

        # Write fragment profile
        frag_profile = read_h5_intervalframe(f["frag_profile"])

        # Read l nfr enrichment values
        l_nfr_enrichment = {}
        nfr_l = f["l_nfr_enrichment"]
        for title in list(nfr_l.keys()):
            title_group = nfr_l[title]
            l_nfr_enrichment[title] = {}
            for chrom in list(title_group.keys()):
                title_group_chrom = title_group[chrom]
                l_nfr_enrichment[title][chrom] = pd.Series(np.array(title_group_chrom["nfr"]),
                                                           index=np.array(title_group_chrom["index"]).astype(str))
                
