import os
import glob
import numpy as np
import pandas as pd
from ailist import AIList
import h5py


def printProgressBar(iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', print_end = "\r"):
    """
    Call in a loop to create terminal progress bar
    
    Parameters
    ----------
        iteration
            int (current iteration)
        total
            int (total iterations)
        prefix
            str (prefix string [defualt:''])
        suffix
            str (suffix string [defualt:''])
        decimals
            int (positive number of decimals in percent complete [defualt:1])
        length
            int (character length of bar [defualt:100])
        fill
            str (bar fill character [defualt:'█'])
        print_end
            str (end character [defualt:'\r'])

    Returns
    -------
        None
    """
    
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print('\r%s |%s| %s%% %s' % (prefix, bar, percent, suffix), end = print_end)
    
    # Print New Line on Complete
    if iteration == total: 
        print()


def find_h5_filenames(names_list, directory="./"):
    """
    Search directory for h5 file names

    Parameters
    ----------
        names_list
            list
        directory
            str

    Returns
    -------
        h5_filenames
            list
    """

    # Record files in directory
    directory_files = glob.glob(os.path.join(directory, "*.h5"), recursive=True)

    # Iterate over names in list
    h5_filenames = []
    for name in names_list:
        found = False
        for dir_fn in directory_files:
            if name in dir_fn:
                h5_filenames.append(dir_fn)
                found = True
                break
        # Could not find file
        if found == False:
            raise FileNotFoundError("Could not find file for " + name)
    
    return h5_filenames


######################################
#     Write methods
######################################


def convert_recarray(rec_array):
    """
    Convert the dtypes of np.recarray to h5py compatible dtypes

    Parameters
    ----------
        rec_array
            numpy.recarray (Array to convert dtypes for)

    Returns
    -------
        rec_array
            numpy.recarray (Array with h5py compatible dtypes)
    """

    dtypes = {o:rec_array.dtype[o].str for o in rec_array.dtype.fields}
    for key_dtype in dtypes:
        if dtypes[key_dtype] == "|O":
            lengths = []
            for i in range(len(rec_array[key_dtype])):
                if isinstance(rec_array[key_dtype][i], str):
                    lengths.append(len(rec_array[key_dtype][i]))
                elif pd.isnull(rec_array[key_dtype][i]):
                    continue
                elif isinstance(rec_array[key_dtype][i], int):
                    continue
                elif isinstance(rec_array[key_dtype][i], float):
                    continue
                else:
                    continue
            max_length = max(lengths) if len(lengths) > 0 else 0
            if max_length > 0:
                dtypes[key_dtype] = "<S{}".format(max_length)
            else:
                dtypes[key_dtype] = "<f8"

    rec_array = rec_array.astype(list(dtypes.items()))

    return rec_array


def write_h5_DataFrame(h5_group, df, axis=0):
    """
    Write pandas.DataFrame to h5py group

    Parameters
    ----------
        h5_group
            h5py.group
        df
            pandas.DataFrame
        axis
            int

    Returns
    -------
        None
    """
    
    # Record axis
    h5_group["axis"] = axis
    h5_group["shape"] = np.array(df.shape)
    
    # Iterate over columns
    if axis == 1:
        if df.index.dtype.kind == "i":
            index_dtypes = "i"
        else:
            index_dtypes = "<S{}".format(df.index.str.len().max())
        if df.columns.dtype.kind == "i":
            columns_dtypes = "<S{}".format(df.columns.str.len().max())
        else:
            columns_dtypes = "<S{}".format(df.columns.str.len().max())

        rec_array = df.T.to_records(index_dtypes=index_dtypes)

        rec_array = convert_recarray(rec_array)

        h5_group["values"] = rec_array
    
    # Iterate over index
    if axis == 0:
        if df.index.dtype.kind == "i":
            index_dtypes = "i"
        else:
            index_dtypes = "<S{}".format(df.index.str.len().max())
        if df.columns.dtype.kind == "i":
            columns_dtypes = "i"
        else:
            columns_dtypes = "<S{}".format(df.columns.str.len().max())

        rec_array = df.to_records(index_dtypes=index_dtypes)

        rec_array = convert_recarray(rec_array)

        h5_group["values"] = rec_array


def write_h5_list_df(h5_group, list_df, axis=0):
    """
    Write list of pd.DataFrames to h5 file

    Parameters
    ----------
        h5_group
            h5py.group
        list_df
            list
        axis
            int

    Returns
    -------
        None
    """

    # Iterate over list
    for i in range(len(list_df)):
        h5_index_group = h5_group.create_group("index_"+str(i), track_order=True)
        write_h5_DataFrame(h5_index_group, list_df[i], axis=axis)


def write_summary(summary, output_h5):
    """
    Write cfDNA_summary object

    Parameters
    ----------
        summary
            dict
        output_h5
            str
    
    Returns
    -------
        None
    """

    # Open h5 file
    output = h5py.File(output_h5, "w")

    # Write descriptions
    output["sam"] = np.array(summary["sam"]).astype(bytes)
    output["n_frags"] = np.array(summary["n_frags"])
    output["bin_size"] = np.array(summary["bin_size"])
    output["genome"] = np.array(summary["genome"]).astype(bytes)
    output["length_dist"] = summary["length_dist"]
    output["tss_profile"] = summary["tss_profile"]
    output["nfr"] = np.array(summary["nfr"])
    output["median_var"] = np.array(summary["median_var"])

    # Write fragments profile
    h5_frag_profile = output.create_group("frag_profile", track_order=True)
    write_h5_DataFrame(h5_frag_profile, summary["frag_profile"], axis=0)

    # Write CNV loglik
    h5_loglik = output.create_group("hmm_loglik", track_order=True)
    write_h5_DataFrame(h5_loglik, summary["hmm_loglik"], axis=0)

    # Write CNV segments
    h5_segments = output.create_group("segments_df", track_order=True)
    write_h5_DataFrame(h5_segments, summary["segments_df"], axis=0)

    # Write CNV bins
    h5_bins = output.create_group("cnv_bins", track_order=True)
    write_h5_DataFrame(h5_bins, summary["cnv_bins"].loc[:,["seqnames","start","end","ratios"]], axis=0)

    # Write WPS bins
    h5_wps_bins = output.create_group("wps_bins", track_order=True)
    write_h5_DataFrame(h5_wps_bins, summary["wps_bins"], axis=0)

    # Write WPS peak dist
    h5_peak_dist = output.create_group("peak_dist", track_order=True)
    write_h5_DataFrame(h5_peak_dist, summary["peak_dist"], axis=0)

    # Write TSS nfr
    h5_tss_nfr = output.create_group("tss_nfr", track_order=True)
    write_h5_DataFrame(h5_tss_nfr, summary["tss_nfr"].loc[:,["seqnames","start","end","nfr"]], axis=0)

    # Write TSS pos_score
    h5_pos_score = output.create_group("tss_pos_score", track_order=True)
    write_h5_DataFrame(h5_pos_score, summary["tss_nfr"].loc[:,["seqnames","start","end","pos_score"]], axis=0)

    # Write chrom shift
    h5_chrom_shift = output.create_group("chrom_shift", track_order=True)
    chrom_shift = pd.DataFrame(list(summary["chrom_shift"].values()),
                                index=list(summary["chrom_shift"].keys()))
    write_h5_DataFrame(h5_chrom_shift, chrom_shift, axis=0)

    # Write genome_len
    h5_genome_len = output.create_group("genome_len", track_order=True)
    genome_len = pd.DataFrame(list(summary["genome_len"].values()),
                                index=list(summary["genome_len"].keys()))
    write_h5_DataFrame(h5_genome_len, genome_len, axis=0)

    # Write WPS peaks
    h5_peaks = output.create_group("wps_peaks", track_order=True)
    for chrom in summary["wps_peaks"]:
        chrom_peaks = np.zeros((summary["wps_peaks"][chrom].size, 2), dtype=int)
        for i, peak in enumerate(summary["wps_peaks"][chrom]):
            chrom_peaks[i,0] = peak.start
            chrom_peaks[i,1] = peak.end
        h5_peaks[chrom] = chrom_peaks

    # Close file
    output.close()


def write_multi_summary(multi_summary, output_h5):
    """
    Write h5 multi_summary metrics

    Parameters
    ----------
        multi_summary
            dict
        output_h5
            str

    Returns
    -------
        None
    """

    # Open h5 file
    output = h5py.File(output_h5, "w")

    # Write descriptions
    output["sam"] = np.array(multi_summary["sam"]).astype(bytes)
    output["n_frags"] = np.array(multi_summary["n_frags"])
    output["bin_size"] = np.array(multi_summary["bin_size"])
    output["genome"] = np.array(multi_summary["genome"]).astype(bytes)
    output["length_dist"] = multi_summary["length_dist"]
    output["tss_profile"] = multi_summary["tss_profile"]
    output["nfr"] = np.array(multi_summary["nfr"])
    output["median_var"] = np.array(multi_summary["median_var"])

    # Write fragments profile
    h5_frag_profile = output.create_group("frag_profile", track_order=True)
    write_h5_DataFrame(h5_frag_profile, multi_summary["frag_profile"], axis=0)

    # Write CNV bins
    h5_bins = output.create_group("cnv_bins", track_order=True)
    write_h5_DataFrame(h5_bins, multi_summary["cnv_bins"], axis=0)

    # Write WPS bins
    h5_wps_bins = output.create_group("wps_bins", track_order=True)
    write_h5_DataFrame(h5_wps_bins, multi_summary["wps_bins"], axis=0)

    # Write WPS peak dist
    h5_peak_dist = output.create_group("peak_dist", track_order=True)
    write_h5_DataFrame(h5_peak_dist, multi_summary["peak_dist"], axis=0)

    # Write TSS nfr
    h5_tss_nfr = output.create_group("tss_nfr", track_order=True)
    write_h5_DataFrame(h5_tss_nfr, multi_summary["tss_nfr"], axis=0)

    # Write TSS pos_score
    h5_pos_score = output.create_group("tss_pos_score", track_order=True)
    write_h5_DataFrame(h5_pos_score, multi_summary["tss_nfr"], axis=0)

    # Write chrom shift
    h5_chrom_shift = output.create_group("chrom_shift", track_order=True)
    chrom_shift = pd.DataFrame(list(multi_summary["chrom_shift"].values()),
                                index=list(multi_summary["chrom_shift"].keys()))
    write_h5_DataFrame(h5_chrom_shift, chrom_shift, axis=0)

    # Write CNV loglik
    h5_loglik = output.create_group("hmm_loglik", track_order=True)
    write_h5_list_df(h5_loglik, multi_summary["hmm_loglik"], axis=0)

    # Write CNV segments
    h5_segments = output.create_group("segments_df", track_order=True)
    write_h5_list_df(h5_segments, multi_summary["segments_df"], axis=0)

    # Close file
    output.close()


######################################
#     Read methods
######################################


def read_h5_DataFrame(h5_group):
    """
    Read pandas.DataFrame from h5py group

    Parameters
    ----------
        h5_group
            h5py.Group

    Returns
    -------
        df
            pandas.DataFrame
    """

    # Define dtype categories
    numeric_dtypes = set(['i','u','b','c'])
    string_dtypes = set(['O','U','S'])
    
    # Record axis
    axis = int(np.array(h5_group["axis"]))
    shape = np.array(h5_group["shape"])
    
    # Read index
    if axis == 0:
        df = pd.DataFrame.from_records(np.array(h5_group["values"]), index="index")
        if df.index.dtype.kind in string_dtypes:
            df.index = df.index.values.astype(str)
        
    # Read columns
    elif axis == 1:
        df = pd.DataFrame.from_records(np.array(h5_group["values"]), index="index").T
        if df.columns.dtype.kind in string_dtypes:
            df.columns = df.columns.values.astype(str)

    # Convert numpy object to str
    for i, dtype in enumerate(df.dtypes):
        if dtype.kind == "O":
            df.iloc[:,i] = df.iloc[:,i].values.astype(str)
        
    return df


def read_h5_list_df(h5_group):
    """
    Read list of pd.DataFrames to h5 file

    Parameters
    ----------
        h5_group
            h5py.Group

    Returns
    -------
        record
            list
    """

    # Iterate over list
    record = []
    for key in list(h5_group.keys()):
        df = read_h5_DataFrame(h5_group[key])
        record.append(df)

    return record


def read_h5_dict_df(h5_group):
    """
    Read dictionary of pd.DataFrames to h5 file

    Parameters
    ----------
        h5_group
            h5py.Group

    Returns
    -------
        record
            list
    """

    # Iterate over list
    record = {}
    for key in list(h5_group.keys()):
        df = read_h5_DataFrame(h5_group[key])
        record[h5_group[key].name] = df

    return record


def read_summary(h5_filename):
    """
    Read summary formatted h5 file to dictionary

    Parameters
    ----------
        h5_filename
            str

    Returns
    -------
        metrics
            dict
    """

    # Open h5 file
    f = h5py.File(h5_filename, "r")

    # Initialize metrics
    metrics = {}
    metrics["sam"] = str(np.array(f["sam"]).astype(str))
    metrics["n_frags"] = int(np.array(f["n_frags"]))
    metrics["bin_size"] = int(np.array(f["bin_size"]))
    metrics["median_var"] = float(np.array(f["median_var"]))
    metrics["genome"] = str(np.array(f["genome"]).astype(str))
    metrics["length_dist"] = np.array(f["length_dist"])
    metrics["nfr"] = tuple(np.array(f["nfr"]))
    metrics["tss_profile"] = np.array(f["tss_profile"])

    # Read DataFrames
    metrics["frag_profile"] = read_h5_DataFrame(f["frag_profile"])
    metrics["hmm_loglik"] = read_h5_DataFrame(f["hmm_loglik"])
    metrics["segments_df"] = read_h5_DataFrame(f["segments_df"])
    metrics["peak_dist"] = read_h5_DataFrame(f["peak_dist"])
    metrics["tss_nfr"] = read_h5_DataFrame(f["tss_nfr"])
    metrics["tss_pos_score"] = read_h5_DataFrame(f["tss_pos_score"])
    metrics["wps_bins"] = read_h5_DataFrame(f["wps_bins"])
    metrics["cnv_bins"] = read_h5_DataFrame(f["cnv_bins"])
    # Remove nans
    metrics["cnv_bins"] = metrics["cnv_bins"].loc[~pd.isnull(metrics["cnv_bins"].loc[:,"ratios"].values),:]
    
    # Read DataFrame to dict
    metrics["chrom_shift"] = read_h5_DataFrame(f["chrom_shift"]).to_dict()["0"]
    metrics["genome_len"] = read_h5_DataFrame(f["genome_len"]).to_dict()["0"]

    # Read Intervals
    metrics["wps_peaks"] = {}
    for chrom in list(f["wps_peaks"].keys()):
        metrics["wps_peaks"][chrom] = AIList()
        intervals = np.array(f["wps_peaks"][chrom])
        metrics["wps_peaks"][chrom].from_array(intervals[0,:], intervals[1,:],
                                               np.arange(intervals.shape[0]),
                                               np.zeros(intervals.shape[0]))
        
    # Close h5 file
    f.close()

    return metrics


def append_metrics(metrics1, metrics2):
    """
    Append two cfDNA_summary h5 file dictionaries

    Parameters
    ----------
        metrics1
            dict
        metrics2
            dict

    Returns
    -------
        total_metrics
            dict
    """

    # Initialize total_metrics
    total_metrics = {}

    # Find sam filenames
    if isinstance(metrics1["sam"], list):
        sam_fn1 = metrics1["sam"]
    else:
        sam_fn1 = os.path.split(metrics1["sam"])[-1]
    sam_fn2 = os.path.split(metrics2["sam"])[-1]

    # Merge bin_size
    if metrics1["bin_size"] != metrics2["bin_size"]:
        raise AttributeError("Bin sizes are not the same")
    else:
        total_metrics["bin_size"] = metrics1["bin_size"]

    # Merge genome
    if metrics1["genome"] != metrics2["genome"]:
        raise AttributeError("Genomes used are not the same")
    else:
        total_metrics["genome"] = metrics1["genome"]

    # Merge attributes
    for key in ["sam", "n_frags", "nfr", "median_var"]:
        if isinstance(metrics1[key], list):
            total_metrics[key] = metrics1[key].copy()
        else:
            total_metrics[key] = [metrics1[key]]
        total_metrics[key].append(metrics2[key])

    # Merge distributions
    for key in ["length_dist", "tss_profile"]:
        if metrics1[key].ndim != 2:
            total_metrics[key] = metrics1[key][np.newaxis]
        else:
            total_metrics[key] = metrics1[key]
        length = min(total_metrics[key].shape[1], len(metrics2[key]))
        total_metrics[key] = np.append(total_metrics[key][:,:length],
                                        metrics2[key][:length][np.newaxis],
                                        axis=0)

    # Merge DataFrames
    for key in ["frag_profile", "peak_dist", "tss_nfr", "tss_pos_score", "wps_bins", "cnv_bins"]:
        # Change column names if first DataFrame
        if metrics1[key].shape[1] == 4:
            metrics1[key].columns = ["seqnames", "start", "end", sam_fn1]
        # Take first two columns
        m1 = metrics1[key].loc[:,["seqnames","start"]].values
        m2 = metrics2[key].loc[:,["seqnames","start"]].values
        # Append first two columns with _
        a = np.array([m1[i,0]+"_"+str(m1[i,1]) for i in range(m1.shape[0])])
        b = np.array([m2[i,0]+"_"+str(m2[i,1]) for i in range(m2.shape[0])])
        # Find which are in common
        in_ab = np.in1d(a, b)
        in_ba = np.in1d(b, a)
        # Set metrics values
        total_metrics[key] = metrics1[key].loc[in_ab,:].copy()
        total_metrics[key].loc[:,sam_fn2] = metrics2[key].loc[in_ba,:].values[:,-1]

    # Append hmm_loglik
    for key in ["hmm_loglik", "segments_df"]:
        if isinstance(metrics1[key], list):
            total_metrics[key] = metrics1[key].copy()
        else:
            total_metrics[key] = [metrics1[key]]
        total_metrics[key].append(metrics2[key])

    return total_metrics


def read_multi_summary(data, directory=".", verbose=False):
    """
    Read multi_summary h5 file

    Parameters
    ----------
        data
            list or str
        directory
            str
        verbose
            bool

    Returns
    -------
        metrics
            dict
    """

    # List provided
    if isinstance(data, list):
        # Determine h5_filenames
        h5_filenames = find_h5_filenames(data, directory)

        # Initialize progress bar
        if verbose:
            progress_length = len(h5_filenames)
            printProgressBar(0, progress_length, prefix='Progress:', suffix='Complete', length=50)

        # Read first h5 file
        metrics = read_summary(h5_filenames[0])
        # Store chrom_shift, only want 1
        chrom_shift = metrics["chrom_shift"]

        # Update progress bar
        if verbose:
            printProgressBar(1, progress_length, prefix='Progress:', suffix='Complete', length=50)

        # Iterate over other h5 files
        for i, h5_filename in enumerate(h5_filenames[1:]):
            metrics_tmp = read_summary(h5_filename)

            # Merge metrics
            metrics = append_metrics(metrics, metrics_tmp)

            # Update Progress Bar
            if verbose:
                printProgressBar(i + 2, progress_length, prefix='Progress:', suffix='Complete', length=50)

        # Re-assign chrom_shift
        metrics["chrom_shift"] = chrom_shift

    # h5 file name provided
    elif isinstance(data, str):

        # Open h5 file
        f = h5py.File(data, "r")

        # Initialize metrics
        metrics = {}
        metrics["sam"] = np.array(f["sam"]).astype(str)
        metrics["n_frags"] = np.array(f["n_frags"])
        metrics["bin_size"] = np.array(f["bin_size"])
        metrics["median_var"] = np.array(f["median_var"])
        metrics["genome"] = np.array(f["genome"]).astype(str)
        metrics["length_dist"] = np.array(f["length_dist"])
        metrics["nfr"] = np.array(f["nfr"])
        metrics["tss_profile"] = np.array(f["tss_profile"])

        # Read DataFrames
        metrics["frag_profile"] = read_h5_DataFrame(f["frag_profile"])
        metrics["peak_dist"] = read_h5_DataFrame(f["peak_dist"])
        metrics["tss_nfr"] = read_h5_DataFrame(f["tss_nfr"])
        metrics["tss_pos_score"] = read_h5_DataFrame(f["tss_pos_score"])
        metrics["wps_bins"] = read_h5_DataFrame(f["wps_bins"])
        metrics["cnv_bins"] = read_h5_DataFrame(f["cnv_bins"])
        
        # Read DataFrame to dict
        metrics["chrom_shift"] = read_h5_DataFrame(f["chrom_shift"]).to_dict()["0"]

        # Read list of DataFrames
        metrics["hmm_loglik"] = read_h5_list_df(f["hmm_loglik"])
        metrics["segments_df"] = read_h5_list_df(f["segments_df"])
            
        # Close h5 file
        f.close()

    return metrics

