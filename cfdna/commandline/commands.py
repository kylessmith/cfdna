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
        cfdna_object = cfDNA(genome_version=args.genome)
        
        # Run CNV for hmm models
        cfdna_object = ngs.segment.cnv_pipeline.call_cnv_pipeline(cfdna_object,
                                                        frags,
                                                        genome_version=args.genome,
                                                        cnv_binsize=args.bin_size,
                                                        hmm_binsize=args.hmm_bin_size,
                                                        nthreads = nthreads)
        
        # Write seg file
        if args.segs:
            seg_fn = prefix + ".segs"
            sample = os.path.split(prefix)[-1]
            ngs.segment.cnv_utilities.write_seg_file(cfdna_object,
                                                        seg_fn,
                                                        sample)

        # Plot
        ngs.plot.cnv_summary(cfdna_object,
                            list(cfdna_object.obs)[0],
                            show = False,
                            save = prefix+"_cnv_plot.png")
