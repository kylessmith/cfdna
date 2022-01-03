import os
import glob
import numpy as np
import pandas as pd
import h5py
from ..core.cfDNA import cfDNA
from ..io.read.read_h5 import read_h5
from ..io.write.write_h5 import write_h5


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


def update_h5(reference_h5, query_h5):
    """
    """

    # Open files
    ref = h5py.File(reference_h5, "r+")
    query = h5py.File(query_h5, "r")

    


def copy_h5(h5_filename, output_h5):
    """
    """

    # Create new h5 file to copy to
    #out_f = h5py.File(output_h5, "w")
    #original_f = h5py.File(h5_filename)

    # Read data
    orig_cfdna = cfDNA()
    read_h5(orig_cfdna, h5_filename)

    # Write data
    write_h5(orig_cfdna, output_h5)


def merge_cfdna_objects(cfdna_object1, cfdna_object2):
    """
    """

    # Merge anno
    cfdna_object1._anno = pd.concat([cfdna_object1._anno, cfdna_object2._anno]).drop_duplicates()

    # Merge intervals
    for key in cfdna_object2.intervals:
        cfdna_object1.add_intervals(key, cfdna_object2.intervals[key].columns.values, cfdna_object2.intervals[key])

    # Merge obs intervals
    for key in cfdna_object2.obs_intervals:
        for obs in cfdna_object2.obs_intervals[key]:
            cfdna_object1.add_obs_intervals(key, obs, cfdna_object2.obs_intervals[key][obs])
    
    # Merge obs values
    for key in cfdna_object2.obs_values:
        cfdna_object1.add_obs_values(key, cfdna_object2.obs_values[key].index.values, cfdna_object2.obs_values[key])

    # Merge HMM states
    for obs in cfdna_object2.hmm_states:
        cfdna_object1.hmm_states[obs] = cfdna_object2.hmm_states[obs]

    # Merge filenames
    cfdna_object1.filenames.update(cfdna_object2.filenames)


def merge_h5(h5_filenames, output_h5, progress_bar=False):
    """
    """

    # Copy first into output
    copy_h5(h5_filenames[0], output_h5)

    # Iterate over h5 files
    for i, h5_filename in enumerate(h5_filenames[1:]):
        # Progress bar
        if progress_bar:
            printProgressBar(i+1, len(h5_filenames), "h5 merge")

        # Read
        cfdna_object1 = cfDNA()
        read_h5(cfdna_object1, output_h5)
        cfdna_object2 = cfDNA()
        read_h5(cfdna_object2, h5_filename)
        # Merge
        merge_cfdna_objects(cfdna_object1, cfdna_object2)
        # Write
        write_h5(cfdna_object1, output_h5)