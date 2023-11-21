import pandas as pd
import numpy as np
import h5py
from projectframe import ProjectFrame


def read_h5(h5_filename: str) -> ProjectFrame:
    """
    Read project frame from HDF5 file

    Parameters
    ----------
        h5_filename : str
            Name of HDF5 file to read from

    Returns
    -------
        pf : :class: `project_frame.ProjectFrame`
            Project frame
    """

    # Read project frame from HDF5 file
    pf = ProjectFrame.read_h5(h5_filename)

    return pf