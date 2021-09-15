from ..core import cfDNA
from .. import Plot


def cfDNA_summarize(args):
    """
    """

    # Detect if multiple bams
    if isinstance(args.bam, list):
        for bam in args.bam:
            # Define prefix
            if args.prefix == "":
                prefix = os.path.split(bam)[-1]
            else:
                prefix = args.prefix + os.path.split(bam)[-1]

            # Create cfDNA object
            cfDNA_object = cfDNA(sam_fn=bam, proportion=args.proportion,
                                 verbose=True, n_jobs=args.n_jobs,
                                 min_size=args.min_length, max_size=args.max_length)

            # Check downsampling
            if args.n_frags != 0:
                cfDNA_object = cfDNA_object.downsample(n_frags=args.n_frags)
            elif args.proportion != 1.0:
                cfDNA_object = cfDNA_object.downsample(args.proportion)

            # Run CNV for hmm models
            cnv = cfDNA_object.call_cnvs(chrom=None, bin_size=1000000, min_length=None, max_length=None,
                                         genome="hg19", bin_bias_h5_fn=None, genome_fn=None, centro_fn=None,
                                         detail_fn=None, outlier_smooth=True, gauss_smooth=False,
                                         normal=[0.1, 0.5, 0.9], ploidy=[1, 2, 3], estimatePloidy=True,
                                         minSegmentBins=25, maxCN=7, bcp_cutoff=0.3)

            # Summarize metrics
            metrics = cfDNA_object.summarize(bin_size=args.bin_size, genome=args.genome, n_threads=args.n_jobs)

            # Relpace HMM metrics
            metrics.metrics["hmm_loglik"] = cnv._segments.hmm_loglik

            # Plot metrics
            print("Saving summary plot:", prefix+"_summary.pdf")
            Plot.summary(metrics, show=False, save=prefix+"_summary.pdf")

            # Write metrics
            print("Saving summary h5:", prefix+"_summary.h5")
            cfDNA_object.write_summary(metrics, prefix+"_summary.h5")
            
            # Write seg file
            if args.segs:
                print("Saving segmentation file:", prefix+".segs")
                seg_fn = prefix + ".segs"
                sample = os.path.split(prefix)[-1]
                metrics.write_seg_file(seg_fn, sample)

    else:
        # Define prefix
        if args.prefix == "":
            prefix = os.path.split(args.bam)[-1]
        else:
            prefix = args.prefix

        # Create cfDNA object
        cfDNA_object = cfDNA(sam_fn=args.bam, bin_size=args.bin_size, proportion=args.proportion,
                            n_frags=args.n_frags, verbose=True, n_jobs=args.n_jobs)

        # Check downsampling
        if args.n_frags != 0:
            cfDNA_object = cfDNA_object.downsample(n_frags=args.n_frags)
        elif args.proportion != 1.0:
            cfDNA_object = cfDNA_object.downsample(args.proportion)

        # Summarize metrics
        metrics = cfDNA_object.summarize(bin_size=args.bin_size, genome=args.genome, n_threads=args.n_jobs)

        # Plot metrics
        Plot.summary(metrics, show=False, save=prefix+"_summary.pdf")
        
        # Write metrics
        cfDNA_object.write_summary(metrics, prefix+"_summary.h5")
        
        # Write seg file
        if args.segs:
            seg_fn = prefix + ".segs"
            sample = os.path.split(prefix)[-1]
            metrics.write_seg_file(seg_fn, sample)


def cfDNA_merge(args):
    """
    Merge multiple summary h5 file into one
    """

    # Detect directory
    directory_split = os.path.split(args.inputs[0])
    if directory_split[0] == "":
        directory = "."
    else:
        directory = directory_split[0]

    # Read multi_summary
    summary = cfDNA_multi_summary(args.inputs, directory=directory, verbose=True)

    # Write multi_summary
    summary.write_h5(args.output)