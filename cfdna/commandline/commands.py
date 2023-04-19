from ..core.core import cfDNA
from ..plot.plot_plt import summary, plot_cnv
from ..io.write.write_text import write_cnv_seg
from ..io.read.read_bam import readBam
from ..processing.cnv.segmentation import call_cnv_pipeline
from ..processing.cnv.segmentation_lite import call_cnv_pipline_lite
from ..processing.summarize.summarize import summarize
from ..utilities.h5_utilities import merge_h5
from ..io.write.write_h5 import write_h5
import os
import multiprocessing
import matplotlib.pyplot as plt


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
        for bam in args.bam:
            if args.prefix == "":
                prefix = os.path.split(bam)[-1]
            else:
                prefix = args.prefix
            # Read bam
            frags = readBam(bam, min_size=args.min_length, max_size=args.max_length,
                                paired=args.single, qcfail=args.noqcfail, mapq_cutoff=args.mapq,
                                proportion=args.proportion, nthreads=nthreads, verbose=args.verbose)

            # Check downsampling
            if args.n_frags != 0:
                frags = frags.downsample(n_frags=args.n_frags)
            elif args.proportion != 1.0:
                frags = frags.downsample(args.proportion)

            # Create cfDNA object
            cfdna_object = cfDNA()
            #cfDNA_object = cfDNA(frags, ref_genome=args.genome, verbose=True)
            
            # Run CNV for hmm models
            call_cnv_pipeline(cfdna_object, frags, cnv_binsize=args.bin_size, hmm_binsize=args.bin_size,
                                          method="bcp_online_both", outlier_smooth=True,
                                          gauss_smooth=False, bcp_cutoff=0.3, normal=[0.1, 0.5, 0.9],
                                          ploidy=[2], estimatePloidy=False, minSegmentBins=25,
                                          maxCN=7, verbose=args.verbose)
            
            # Write seg file
            if args.segs:
                print("Saving segmentation file:", prefix+".segs")
                seg_fn = prefix + ".segs"
                sample = os.path.split(prefix)[-1]
                write_cnv_seg(cfdna_object, sample, seg_fn)

            # Plot
            plot_cnv(cfdna_object, list(cfdna_object.filenames)[0], title=list(cfdna_object.filenames)[0],
                    show=False, save=prefix+"_cnv_plot.png", ax=None)

    else:
        if args.prefix == "":
            prefix = os.path.split(args.bam)[-1]
        else:
            prefix = args.prefix
        # Read bam
        #call_cnv_pipline_lite(args.bam, cnv_binsize=args.bin_size, hmm_binsize=args.bin_size,
        #                                method="bcp_online_both", outlier_smooth=True,
        #                                gauss_smooth=False, bcp_cutoff=0.3, normal=[0.1, 0.5, 0.9],
        #                                ploidy=[2], estimatePloidy=False, minSegmentBins=25,
        #                                maxCN=7, verbose=False)

        frags = readBam(args.bam, min_size=args.min_length, max_size=args.max_length,
                            paired=args.single, qcfail=args.noqcfail, mapq_cutoff=args.mapq,
                            proportion=args.proportion, nthreads=nthreads, verbose=args.verbose)

        # Check downsampling
        if args.n_frags != 0:
            frags = frags.downsample(n_frags=args.n_frags)
        elif args.proportion != 1.0:
            frags = frags.downsample(args.proportion)
        
        # Create cfDNA object
        cfdna_object = cfDNA()
        #cfDNA_object = cfDNA(fragments, ref_genome=args.genome, verbose=True)

        # Run CNV for hmm models
        call_cnv_pipline(cfdna_object, frags, cnv_binsize=args.bin_size, hmm_binsize=args.bin_size,
                                        method="bcp_online_both", outlier_smooth=True,
                                        gauss_smooth=False, bcp_cutoff=0.3, normal=[0.1, 0.5, 0.9],
                                        ploidy=[2], estimatePloidy=False, minSegmentBins=25,
                                        maxCN=7, verbose=args.verbose)

        # Write seg file
        if args.segs:
            seg_fn = prefix + ".segs"
            sample = os.path.split(prefix)[-1]
            write_cnv_seg(cfdna_object, sample, seg_fn)

        # Plot
        plot_cnv(cfdna_object, list(cfdna_object.filenames)[0], title=list(cfdna_object.filenames)[0],
                 show=False, save=prefix+"_cnv_plot.png", ax=None)


def cfDNA_summarize(args):
    """
    """

    # Check number of cores
    if args.nthreads == -1:
        nthreads = multiprocessing.cpu_count()
    else:
        nthreads = args.nthreads

    # Detect if multiple bams
    if isinstance(args.bam, list):
        for bam in args.bam:
            # Define prefix
            if args.prefix == "":
                prefix = os.path.split(bam)[-1]
            else:
                prefix = args.prefix + os.path.split(bam)[-1]

            # Define sample
            sample = os.path.split(prefix)[-1]

            # Create cfDNA object
            cfdna_object = cfDNA()
            frags = readBam(bam, proportion=args.proportion,
                                 verbose=args.verbose, min_size=args.min_length, max_size=args.max_length,
                                 nthreads=nthreads, is_meth=args.meth)

            # Check downsampling
            if args.n_frags != 0:
                frags = frags.downsample(n_frags=args.n_frags)
            elif args.proportion != 1.0:
                frags = frags.downsample(args.proportion)

            # Run summary
            if args.verbose: print("Summarizing", flush=True)
            summarize(cfdna_object, frags, args.bin_size, progress_bar=args.progressbar, verbose=args.verbose)

            # Plot metrics
            if args.verbose: print("Saving summary plot:", prefix+"_summary.pdf")
            summary(cfdna_object, sample, show=False, save=prefix+"_summary.pdf")

            # Write metrics
            if args.verbose: print("Saving summary h5:", prefix+"_summary.h5")
            write_h5(cfdna_object, prefix+"_summary.h5")
            
            # Write seg file
            if args.segs:
                if args.verbose: print("Saving segmentation file:", prefix+".segs")
                seg_fn = prefix + ".segs"
                sample = os.path.split(prefix)[-1]
                write_cnv_seg(cfdna_object, sample, seg_fn)

    else:
        # Define BAM
        bam = args.bam

        # Define prefix
        if args.prefix == "":
            prefix = os.path.split(bam)[-1]
        else:
            prefix = args.prefix + os.path.split(bam)[-1]

        # Define sample
        sample = os.path.split(prefix)[-1]

        # Create cfDNA object
        cfdna_object = cfDNA()
        frags = readBam(bam, proportion=args.proportion,
                                verbose=True, min_size=args.min_length, max_size=args.max_length,
                                nthreads=nthreads, is_meth=args.meth)

        # Check downsampling
        if args.n_frags != 0:
            frags = frags.downsample(n_frags=args.n_frags)
        elif args.proportion != 1.0:
            frags = frags.downsample(args.proportion)

        # Run summary
        if args.verbose: print("Summarizing", flush=True)
        summarize(cfdna_object, frags, args.bin_size, progress_bar=args.progressbar, verbose=args.verbose)

        # Plot metrics
        if args.verbose: print("Saving summary plot:", prefix+"_summary.pdf")
        summary(cfdna_object, sample, show=False, save=prefix+"_summary.pdf")

        # Write metrics
        if args.verbose: print("Saving summary h5:", prefix+"_summary.h5")
        write_h5(cfdna_object, prefix+"_summary.h5")
        
        # Write seg file
        if args.segs:
            if args.verbose: print("Saving segmentation file:", prefix+".segs")
            seg_fn = prefix + ".segs"
            sample = os.path.split(prefix)[-1]
            write_cnv_seg(cfdna_object, sample, seg_fn)

    # Clean up matplotlib plots
    plt.close("all")


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

    # Merge h5 files
    merge_h5(args.inputs, args.output, progress_bar=True)