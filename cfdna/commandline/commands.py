from ..core.core import cfDNA
import os
import multiprocessing
import matplotlib.pyplot as plt
import ngsfragments as ngs


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
                cfdna_object = ngs.segment.cnv_pipeline.call_cnv_pipeline(cfdna_object,
                                                                frags,
                                                                genome_version=args.genome,
                                                                cnv_binsize=args.bin_size,
                                                                hmm_binsize=args.hmm_bin_size,
                                                                nthreads = nthreads,
                                                                use_normal = args.use_normal,
                                                                keep_sex_chroms = args.add_sex,
                                                                normal = [0.1, 0.25, 0.5, 0.75, 0.9],
                                                                ploidy = [2,3],
                                                                estimatePloidy = True,
                                                                scStates = [1, 3],
                                                                minSegmentBins = 25,
                                                                maxCN = 5)
            else:
                cfdna_object = ngs.segment.cnv_pipeline.call_cnv_pipeline(cfdna_object,
                                                            frags,
                                                            genome_version=args.genome,
                                                            cnv_binsize=args.bin_size,
                                                            hmm_binsize=args.hmm_bin_size,
                                                            nthreads = nthreads,
                                                            use_normal = args.use_normal,
                                                            keep_sex_chroms = args.add_sex,
                                                            normal = [0.1, 0.25, 0.5, 0.75, 0.9],
                                                            ploidy = [2,3],
                                                            estimatePloidy = True,
                                                            scStates = None,
                                                            minSegmentBins = 25,
                                                            maxCN = 5)
            
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
        gene_activity = ngs.metrics.gene_activity(frags,
                                                genome_version=args.genome,
                                                feature=args.feature,
                                                min_length=args.min_length,
                                                max_length=args.max_length,
                                                verbose=args.verbose)
        
        # Correct gene activity
        gene_activity = ngs.metrics.correct_gene_activity(frags,
                                                        gene_activity,
                                                        correct_cnv = False,
                                                        genome_version = args.genome,
                                                        feature = args.feature,
                                                        verbose = args.verbose)

        # Add to cfDNA object
        cfdna_object.add_values("gene_activity", gene_activity)

        # Write h5
        if args.cache != "":
            cfdna_object.to_h5(args.cache)
        else:
            cfdna_object.to_h5(prefix+".h5")

        