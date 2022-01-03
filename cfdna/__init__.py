from __future__ import absolute_import
from .core.cfDNA import cfDNA
from . import plot as pl
from . import io
from . import commandline as cmd
from . import utilities as utils
from . import processing as proc

from ngsfragments import fragments


# This is extracted automatically by the top-level setup.py.
__version__ = '1.0.0'

__author__ = "Kyle S. Smith"

__doc__ = """\

API
======

Central class
-------------

.. autosummary::
   :toctree: .
   
   cfDNA
   
   
Utility functions
-----------------

.. autosummary::
   :toctree: .
   
   utils.h5_utilities


Plotting functions
------------------
    
.. autosummary::
   :toctree: .
   
   pl.plot_plt
   pl.plot_bokeh


IO functions
------------------
    
.. autosummary::
   :toctree: .
   
   io.write
   io.read

"""