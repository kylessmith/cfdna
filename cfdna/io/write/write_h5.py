import os
import glob
import numpy as np
import pandas as pd
from projectframe import ProjectFrame


def write_h5(pf: ProjectFrame, h5_filename: str) -> None:
    """
    Write project frame to HDF5 file

    Parameters
    ----------
        pf : :class: `project_frame.ProjectFrame`
            Project frame
        h5_filename : str
            Name of HDF5 file to write to

    Returns
    -------
        None
    """

    # Write project frame to HDF5 file
    pf.write_h5(h5_filename)

    return None
