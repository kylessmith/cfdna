from ..core.core import cfDNA
import os
import multiprocessing
import matplotlib.pyplot as plt
import numpy as np
import ngsfragments as ngs
from intervalframe import IntervalFrame
from ailist import LabeledIntervalArray
import genome_info
import glob

# Local imports
from ..data.import_data import get_data_file
from ..tools.nucleosome.wps import wps_gene_fft
from ..tools.nucleosome.wps import predict_nucleosomes


def CNV_calling(args):
    """
    """

    # Check number of cores
    if args.nthreads == -1:
        nthreads = multiprocessing.cpu_count()
    else:
        nthreads = args.nthreads

    # Detect if multiple bams
    if isinstance(args.bam, list):
        bams = args.bam
    else:
        bams = [args.bam]

    # Read PON
    if args.pon != "":
        pon = IntervalFrame.read_parquet(args.pon)

    # Read additional blacklist
    if args.genome == "hg38":
        blacklist = IntervalFrame.read_parquet(get_data_file("cfdna_hg38_blacklist.parquet"))
    else:
        blacklist = None

    for bam in bams:
        if args.prefix == "":
            prefix = os.path.split(bam)[-1].split(".bam")[0]
        else:
            prefix = args.prefix

        # Create cfDNA object
        if args.cache != "":
            # Check if cache exists
            if os.path.exists(args.cache):
                cfdna_object = cfDNA.from_h5(args.cache)
        else:
            cfdna_object = cfDNA(genome_version=args.genome)

            # Read bam
            if args.nanopore:
                import pysam
                samfile = pysam.AlignmentFile(bam, "rb")
                ail = LabeledIntervalArray()
                iters = samfile.fetch(until_eof=True)
                for x in iters:
                    if x.is_unmapped == False and x.is_secondary == False and x.is_supplementary == False:
                        length = x.reference_end - x.reference_start
                        if length >= args.min_length and length <= args.max_length:
                            ail.add(int(x.reference_start), int(x.reference_end), x.reference_name)
                frags = ngs.Fragments(genome_version=args.genome)
                frags.frags = ail
                frags.sam_file = bam
            else:
                frags = ngs.io.from_sam(bam,
                                        genome_version=args.genome,
                                        min_size = args.min_length,
                                        max_size = args.max_length,
                                        paired = not args.single,
                                        qcfail = args.qcfail,
                                        mapq_cutoff = args.mapq,
                                        nthreads=nthreads,
                                        verbose=args.verbose)
            
            # Run CNV for hmm models
            if args.clonal:
                cnvs = ngs.segment.CNVcaller(genome_version=args.genome,
                                            cnv_binsize=args.bin_size,
                                            hmm_binsize=args.hmm_bin_size,
                                            use_normal = args.use_normal,
                                            keep_sex_chroms = args.add_sex,
                                            normal = [0.1, 0.25, 0.5, 0.75, 0.9],
                                            ploidy = [2, 3, 4],
                                            estimatePloidy = True,
                                            scStates = [1, 3],
                                            minSegmentBins = 25,
                                            maxCN = 5)
                
            else:
                cnvs = ngs.segment.CNVcaller(genome_version=args.genome,
                                            cnv_binsize=args.bin_size,
                                            hmm_binsize=args.hmm_bin_size,
                                            use_normal = args.use_normal,
                                            keep_sex_chroms = args.add_sex,
                                            normal = [0.1, 0.25, 0.5, 0.75, 0.9],
                                            ploidy = [2, 3, 4],
                                            estimatePloidy = True,
                                            scStates = None,
                                            minSegmentBins = 25,
                                            maxCN = 5)
            # Call CNVs
            sample = os.path.split(prefix)[-1]
            if args.pon != "":
                cnvs.predict_cnvs(frags, normal_data=pon, additional_blacklist=blacklist, additional_blacklist_cutoff=0.3)
            else:
                cnvs.predict_cnvs(frags, additional_blacklist=blacklist, additional_blacklist_cutoff=0.3)
            cfdna_object = cfDNA(pf=cnvs.pf, genome_version=args.genome)
            cfdna_object.log_fragments(frags)
            cfdna_object.params["cnv_binsize"] = args.bin_size
            # Calculate segment variance
            cfdna_object.obs_intervals[sample]["cnv_segments"].annotate(cfdna_object.intervals["cnv_bins"], sample, "var")
            
            #cfdna_object.add_anno("n_fragments", sample, len(frags.frags))
            # Calculate segment variance
            #cfdna_object.obs_intervals[sample]["cnv_segments"].annotate(cfdna_object.intervals["cnv_bins"], sample, "var")
            #length_dist = frags.length_dist()
            #length_dist.name = sample
            #cfdna_object.add_values("length_dist", length_dist)
            
            # Add WPS
            if args.add_wps:
                scores = ngs.metrics.wps_windows(frags,
                                                protection = 120,
                                                min_length = args.min_length,
                                                max_length = args.max_length,
                                                feature = "tss",
                                                smooth = False)
                cfdna_object.add_obs_values(list(cfdna_object.obs)[0], "tss_wps", scores)
        
        # Write seg file
        if args.segs:
            seg_fn = prefix + ".seg"
            sample = os.path.split(prefix)[-1]
            ngs.segment.cnv_utilities.write_seg_file(cfdna_object,
                                                        seg_fn,
                                                        sample)

        # Plot
        ngs.plot.cnv_summary(cfdna_object,
                            list(cfdna_object.obs)[0],
                            add_wps = args.add_wps,
                            show = False,
                            save = prefix+"_cnv_plot.pdf")
        
        # Write h5
        if args.cache != "" and os.path.exists(args.cache) == False:
            cfdna_object.to_h5(args.cache)

        # Write annotations
        if args.anno_file:
            cfdna_object.anno.engine.df.to_csv(prefix+"_metrics.txt", header=True, index=True, sep="\t")
        
        # Write seg annotations
        if args.anno_segs:
            sample = os.path.split(prefix)[-1]
            df = cfdna_object.obs_intervals[sample]["cnv_segments"].df
            df.loc[:,"chrom"] = cfdna_object.obs_intervals[sample]["cnv_segments"].index.labels
            df.loc[:,"start"] = cfdna_object.obs_intervals[sample]["cnv_segments"].index.starts
            df.loc[:,"end"] = cfdna_object.obs_intervals[sample]["cnv_segments"].index.ends
            df.loc[:,"sample"] = sample
            df.loc[:,"n_bins"] = ((df.loc[:,"end"].values - df.loc[:,"start"].values) / args.bin_size).astype(int)
            df = df.loc[:,['sample', 'chrom', 'start', 'end', 'copy_number', 'event', 'subclone_status',
                           'logR_Copy_Number', 'Corrected_Copy_Number', 'Corrected_Call', 'var', 'n_bins', 'median']]
            df.to_csv(prefix+"_seg_annotations.seg", header=True, index=False, sep="\t")

            # Drop columns
            df.drop(columns=["chrom", "start", "end"], inplace=True)

        # Write bins
        if args.bins_file:
            sample = os.path.split(prefix)[-1]
            df = cfdna_object.intervals["cnv_bins"].df
            df.loc[:,"chrom"] = cfdna_object.intervals["cnv_bins"].index.labels
            df.loc[:,"start"] = cfdna_object.intervals["cnv_bins"].index.starts
            df.loc[:,"end"] = cfdna_object.intervals["cnv_bins"].index.ends
            df.loc[:,"ratio"] = df.loc[:,sample].values
            df = df.loc[:,['chrom', 'start', 'end', 'ratio']]
            df.to_csv(prefix+"_bins.txt", header=True, index=False, sep="\t")

            # Drop columns
            df.drop(columns=["chrom", "start", "end", "ratio"], inplace=True)


def GeneActivity(args):
    """
    """

    # Check number of cores
    if args.nthreads == -1:
        nthreads = multiprocessing.cpu_count()
    else:
        nthreads = args.nthreads

    # Detect if multiple bams
    if isinstance(args.bam, list):
        bams = args.bam
    else:
        bams = [args.bam]

    for bam in bams:
        if args.prefix == "":
            prefix = os.path.split(bam)[-1].split(".bam")[0]
        else:
            prefix = args.prefix
        sample = os.path.split(prefix)[-1]
        
        # Read bam
        frags = ngs.io.from_sam(bam,
                                genome_version=args.genome,
                                min_size = args.min_length,
                                max_size = args.max_length,
                                paired = not args.single,
                                qcfail = args.qcfail,
                                mapq_cutoff = args.mapq,
                                nthreads=nthreads,
                                verbose=args.verbose)

        # Create cfDNA object
        if args.cache != "":
            cfdna_object = cfDNA.from_h5(args.cache)
        else:
            cfdna_object = cfDNA(genome_version=args.genome)

        # Get gene activity
        gene_pattern = wps_gene_fft(frags,
                                    genome_version = args.genome,
                                    protection = args.protection,
                                    min_length = args.min_length,
                                    max_length = args.max_length,
                                    freq_range = (120, 280),
                                    scale = True,
                                    feature = args.feature,
                                    verbose = args.verbose)
        #gene_pattern = gene_pattern.loc[gene_pattern.sum(axis=1).values != 0,:]
        gene_activity = gene_pattern.loc[:,193:205].mean(axis=1).to_frame()
        gene_activity.columns = [sample]

        # Add to cfDNA object
        cfdna_object.add_values("gene_activity", gene_activity)

    # Write parquet
    cfdna_object.values["gene_activity"].to_parquet("gene_activity.parquet")


def GenomeCoverage(args):
    """
    """
    import pyBigWig

    # Check number of cores
    if args.nthreads == -1:
        nthreads = multiprocessing.cpu_count()
    else:
        nthreads = args.nthreads

    # Detect if multiple bams
    if isinstance(args.bam, list):
        bams = args.bam
    else:
        bams = [args.bam]

    for bam in bams:
        if args.prefix == "":
            prefix = os.path.split(bam)[-1].split(".bam")[0]
        else:
            prefix = args.prefix
        sample = os.path.split(prefix)[-1]
        
        # Read bam
        frags = ngs.io.from_sam(bam,
                                genome_version=args.genome,
                                min_size = args.min_length,
                                max_size = args.max_length,
                                paired = not args.single,
                                qcfail = args.qcfail,
                                mapq_cutoff = args.mapq,
                                nthreads=nthreads,
                                verbose=args.verbose)
        # Calculate coverage
        bw = pyBigWig.open(sample+"_coverage.bw", "w")
        header = []
        chrom_ranges = frags.frags.label_ranges
        for chrom in chrom_ranges:
            ranges = chrom_ranges[chrom]
            header.append((chrom, ranges[1]))
        bw.addHeader(header)
        chroms = frags.frags.unique_labels
        for chrom in chroms:
            if args.verbose: print(chrom)
            cov = ngs.coverage.coverage(frags, chrom=chrom, min_length=args.min_length, max_length=args.max_length)
            ngs.coverage.normalize_coverage(cov, method="mean")
            cov[chrom] = cov[chrom].astype(np.float32)

            # Write to bigwig
            bw.addEntries(chrom, cov[chrom].index.values, values=cov[chrom].values, span=1)
        
        bw.close()


def GenomeWPS(args):
    """
    """
    import pyBigWig

    # Check number of cores
    if args.nthreads == -1:
        nthreads = multiprocessing.cpu_count()
    else:
        nthreads = args.nthreads

    # Detect if multiple bams
    if isinstance(args.bam, list):
        bams = args.bam
    else:
        bams = [args.bam]

    for bam in bams:
        if args.prefix == "":
            prefix = os.path.split(bam)[-1].split(".bam")[0]
        else:
            prefix = args.prefix
        sample = os.path.split(prefix)[-1]
        
        # Read bam
        frags = ngs.io.from_sam(bam,
                                genome_version=args.genome,
                                min_size = args.min_length,
                                max_size = args.max_length,
                                paired = not args.single,
                                qcfail = args.qcfail,
                                mapq_cutoff = args.mapq,
                                nthreads=nthreads,
                                verbose=args.verbose)
        # Calculate coverage
        bw = pyBigWig.open(sample+"_wps.bw", "w")
        header = []
        chrom_ranges = frags.frags.label_ranges
        for chrom in chrom_ranges:
            ranges = chrom_ranges[chrom]
            header.append((chrom, ranges[1]))
        bw.addHeader(header)
        chroms = frags.frags.unique_labels
        for chrom in chroms:
            if args.verbose: print(chrom)
            cov = ngs.coverage.wps(frags, chrom=chrom, min_length=args.min_length, max_length=args.max_length, protection=args.protection)
            ngs.coverage.normalize_coverage(cov, method="mean")
            cov[chrom] = cov[chrom].astype(np.float32)

            # Write to bigwig
            bw.addEntries(chrom, cov[chrom].index.values, values=cov[chrom].values, span=1)
        
        bw.close()


def WPSpeaks(args):
    """
    """
    import pyBigWig

    # Check number of cores
    if args.nthreads == -1:
        nthreads = multiprocessing.cpu_count()
    else:
        nthreads = args.nthreads

    # Detect if multiple bams
    if isinstance(args.bam, list):
        bams = args.bam
    else:
        bams = [args.bam]

    for bam in bams:
        if args.prefix == "":
            prefix = os.path.split(bam)[-1].split(".bam")[0]
        else:
            prefix = args.prefix
        sample = os.path.split(prefix)[-1]
        
        # Read bam
        frags = ngs.io.from_sam(bam,
                                genome_version=args.genome,
                                min_size = args.min_length,
                                max_size = args.max_length,
                                paired = not args.single,
                                qcfail = args.qcfail,
                                mapq_cutoff = args.mapq,
                                nthreads=nthreads,
                                verbose=args.verbose)
        
        peaks = predict_nucleosomes(frags,
                        protection = args.protection,
                        merge_distance = 5,
				        min_length = args.min_length,
                        max_length = args.max_length)
        iframe_peaks = IntervalFrame(intervals=peaks)
        iframe_peaks.to_bed(prefix+"_nucleosomes.bed")


def create_pon(args):
    """
    """

    # Load genome
    g = genome_info.GenomeInfo(args.genome)

    # Check number of cores
    if args.nthreads == -1:
        nthreads = multiprocessing.cpu_count()
    else:
        nthreads = args.nthreads

    # Detect if multiple bams
    bams = glob.glob(args.dir + "/*.bam")

    # Create pon
    pon = LabeledIntervalArray.create_bin(g["chrom_sizes"], bin_size=args.bin_size)
    pon = IntervalFrame(intervals=pon)
    pon = pon.loc[g["main_chromosomes"],:]

    for bam in bams:
        if args.verbose: print(bam)
        prefix = os.path.split(bam)[-1].split(".bam")[0]
        sample = os.path.split(prefix)[-1]
        
        # Read bam
        frags = ngs.io.from_sam(bam,
                                genome_version=args.genome,
                                min_size = args.min_length,
                                max_size = args.max_length,
                                paired = not args.single,
                                qcfail = args.qcfail,
                                mapq_cutoff = args.mapq,
                                nthreads=nthreads,
                                verbose=args.verbose)

        # Bin coverage
        cnv = ngs.segment.CNVcaller(genome_version=args.genome, cnv_binsize=args.bin_size)
        bins = cnv.bin_data(frags)
        bins = bins.loc[g["main_chromosomes"],:]
        bins = bins.exact_match(pon)
        pon = pon.exact_match(bins)

        # Calculate coverage
        pon.df.loc[:,sample] = bins.df.loc[:,sample].values

    # Write to parquet
    pon.to_parquet(args.out)