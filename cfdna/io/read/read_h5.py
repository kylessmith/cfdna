import pandas as pd
import numpy as np
import h5py
from intervalframe.read.read_h5 import read_h5_intervalframe


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

    # Record columns type
    col_type = np.array(h5_group["columns_type"])[0].decode()

    # Read index
    if axis == 0:
        df = pd.DataFrame.from_records(np.array(h5_group["values"]), index="index")
        if df.index.dtype.kind in string_dtypes:
            df.index = df.index.values.astype(str)
        
    # Read columns
    elif axis == 1:
        if col_type == "single":
            df = pd.DataFrame.from_records(np.array(h5_group["values"]), index="index").T
            if df.columns.dtype.kind in string_dtypes:
                df.columns = df.columns.values.astype(str)
        elif col_type == "multi":
            df = pd.DataFrame.from_records(np.array(h5_group["values"]), index=["level_0","level_1"]).T
            df.columns = pd.MultiIndex.from_tuples([(a.decode(), b.decode()) for a,b in df.columns.values],
                                                    names=["level_0","level_1"])

    # Convert numpy object to str
    for i, dtype in enumerate(df.dtypes):
        if dtype.kind == "O":
            df.iloc[:,i] = df.iloc[:,i].values.astype(str)
        
    return df


def read_h5(cfdna_object, h5_filename):
    """
    """

    # Open file
    f = h5py.File(h5_filename, "r")

    # Read anno
    cfdna_object._anno = read_h5_DataFrame(f["anno"])

    # Read intervals
    for key in list(f["intervals"].keys()):
        cfdna_object._intervals[key] = read_h5_intervalframe(f["intervals"][key])

    # Read obs intervals
    for key in list(f["obs_intervals"]):
        cfdna_object._obs_intervals[key] = {}
        for obs in list(f["obs_intervals"][key]):
            cfdna_object._obs_intervals[key][obs] = read_h5_intervalframe(f["obs_intervals"][key][obs])

    # Read obs values
    for key in list(f["obs_values"]):
        cfdna_object._obs_values[key] = read_h5_DataFrame(f["obs_values"][key])

    # Read chroms
    cfdna_object.chroms = set(np.array(f["chroms"]).astype(str))

    # Read hmm
    hmm_states = {}
    for key in list(f["hmm_states"]):
        hmm_states[key] = {}
        hmm_states[key]["mus"] = np.array(f["hmm_states"][key]["mus"])
        hmm_states[key]["lambdas"] = np.array(f["hmm_states"][key]["lambdas"])
        hmm_states[key]["nu"] = np.array(f["hmm_states"][key]["nu"])[0]
        hmm_states[key]["states"] = pd.Series(np.array(f["hmm_states"][key]["states"]), name="Sample_1")
        hmm_states[key]["subclone"] = pd.Series(np.array(f["hmm_states"][key]["subclone"]), name="Sample_1")
        hmm_states[key]["n"] = np.array(f["hmm_states"][key]["n"])[0]
        hmm_states[key]["phi"] = np.array(f["hmm_states"][key]["phi"])[0]
        hmm_states[key]["sp"] = np.array(f["hmm_states"][key]["sp"])[0]
        hmm_states[key]["Frac_genome_subclonal"] = np.array(f["hmm_states"][key]["Frac_genome_subclonal"])[0]
        hmm_states[key]["cn"] = read_h5_intervalframe(f["hmm_states"][key]["cn"])
    cfdna_object.hmm_states = hmm_states

    # Read file names
    cfdna_object.filenames = set(np.array(f["filenames"]).astype(str))

    # Read bin sizes
    if "cnv_binsize" in list(f.keys()):
        cfdna_object.cnv_binsize = np.array(f["cnv_binsize"])[0]
    if "hmm_binsize" in list(f.keys()):
        cfdna_object.hmm_binsize = np.array(f["hmm_binsize"])[0]
    if "profile_binsize" in list(f.keys()):
        cfdna_object.profile_binsize = np.array(f["profile_binsize"])[0]

    # Read chrom lengths
    chroms = np.array(f["chrom_lengths"]["chroms"]).astype(str)
    lengths = np.array(f["chrom_lengths"]["lengths"])
    cfdna_object.chrom_lengths = {chrom:lengths[i] for i, chrom in enumerate(chroms)}

    # Close file
    f.close()