import os
import glob
import numpy as np
import pandas as pd
import h5py
from intervalframe.write.write_h5 import write_h5_intervalframe


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

    # Handle MultiIndex
    if isinstance(df.columns, pd.MultiIndex):
        h5_group["columns_type"] = np.array([b"multi"])
    else:
        h5_group["columns_type"] = np.array([b"single"])

    # Iterate over columns
    if axis == 1:
        if df.index.dtype.kind == "i":
            index_dtypes = "i"
        else:
            index_dtypes = "<S{}".format(df.index.str.len().max())
        
        #if df.columns.dtype.kind == "i":
        #    columns_dtypes = "<S{}".format(df.columns.str.len().max())
        #else:
        #    columns_dtypes = "<S{}".format(df.columns.str.len().max())

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


def write_h5(cfdna_object, h5_filename):
    """
    """
    
    # Open file
    f = h5py.File(h5_filename, "w")

    # Write anno
    anno = f.create_group("anno")
    write_h5_DataFrame(anno, cfdna_object.anno, axis=1)

    # Write intervals
    intervals = f.create_group("intervals")
    for key in cfdna_object.intervals:
        key_intervals = intervals.create_group(key)
        write_h5_intervalframe(cfdna_object.intervals[key], key_intervals)

    # Write obs intervals
    obs_intervals = f.create_group("obs_intervals")
    for key in cfdna_object.obs_intervals:
        #key_obs_intervals = obs_intervals.create_group(key+"_obs_intervals")
        key_obs_intervals = obs_intervals.create_group(key)
        for obs in cfdna_object.obs_intervals[key]:
            #key_obs_obs_intervals = key_obs_intervals.create_group(key+"_obs_obs_intervals")
            key_obs_obs_intervals = key_obs_intervals.create_group(obs)
            write_h5_intervalframe(cfdna_object.obs_intervals[key][obs], key_obs_obs_intervals)

    # Write obs values
    obs_values = f.create_group("obs_values")
    for key in cfdna_object.obs_values:
        #key_obs_values = obs_values.create_group(key+"obs_values")
        key_obs_values = obs_values.create_group(key)
        write_h5_DataFrame(key_obs_values, cfdna_object.obs_values[key], axis=1)

    # Write chroms
    f["chroms"] = np.array(list(cfdna_object.chroms)).astype(bytes)

    # Write hmm
    hmm_states = f.create_group("hmm_states")
    for key in cfdna_object.hmm_states:
        #key_hmm_states = hmm_states.create_group(key+"_hmm_states")
        key_hmm_states = hmm_states.create_group(key)
        key_hmm_states["mus"] = cfdna_object.hmm_states[key]["mus"]
        key_hmm_states["lambdas"] = cfdna_object.hmm_states[key]["lambdas"]
        key_hmm_states["nu"] = np.array([cfdna_object.hmm_states[key]["nu"]])
        key_hmm_states["states"] = cfdna_object.hmm_states[key]["states"].values
        key_hmm_states["subclone"] = cfdna_object.hmm_states[key]["subclone"].values
        key_hmm_states["n"] = np.array([cfdna_object.hmm_states[key]["n"]])
        key_hmm_states["phi"] = np.array([cfdna_object.hmm_states[key]["phi"]])
        key_hmm_states["sp"] = np.array([cfdna_object.hmm_states[key]["sp"]])
        key_hmm_states["Frac_genome_subclonal"] = np.array([cfdna_object.hmm_states[key]["Frac_genome_subclonal"]])
        #cn_key_hmm_states = hmm_states.create_group("cn_"+key+"_hmm_states")
        cn_key_hmm_states = key_hmm_states.create_group("cn")
        write_h5_intervalframe(cfdna_object.hmm_states[key]["cn"], cn_key_hmm_states)

    # Write file names
    f["filenames"] = np.array(list(cfdna_object.filenames)).astype(bytes)

    # Write bin sizes
    if cfdna_object.cnv_binsize is not None:
        f["cnv_binsize"] = np.array([cfdna_object.cnv_binsize])
    if cfdna_object.hmm_binsize is not None:
        f["hmm_binsize"] = np.array([cfdna_object.hmm_binsize])
    if cfdna_object.profile_binsize is not None:
        f["profile_binsize"] = np.array([cfdna_object.profile_binsize])

    # Write chrom lengths
    chrom_lengths = f.create_group("chrom_lengths")
    chrom_lengths["chroms"] = np.array(list(cfdna_object.chrom_lengths.keys())).astype(bytes)
    chrom_lengths["lengths"] = np.array(list(cfdna_object.chrom_lengths.values()))

    # Close file
    f.close()