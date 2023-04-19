import numpy as np
import pandas as pd


def write_seg(iframe, filename, sample_name, value_key, binsize=100000):
    """
    """
    
    # Create pandas dataframe
    chroms = iframe.index.extract_labels()
    starts = iframe.starts()
    ends = iframe.ends()
    nbins = (ends - starts) / binsize
    seg_csv = pd.DataFrame([np.repeat(sample_name, iframe.shape[0]), chroms,
                            starts, ends, nbins.astype(int), iframe.df.loc[:,value_key]],
                           index=["sample","chrom","start","end","n_bins","value"]).T
    
    # Write
    seg_csv.to_csv(filename, header=True, index=None, sep="\t")


def write_cnv_seg(cfdna_object, key, filename):
    """
    """

    # Write seg file
    write_seg(cfdna_object.obs_intervals["cnv_segments"][key],
              filename, key, "median", cfdna_object.cnv_binsize)