from __future__ import absolute_import
from .core.core import cfDNA
from . import plot as pl
from . import io
from . import commandline as cmd
from . import tools as tl

from ngsfragments import Fragments


# This is extracted automatically by the top-level setup.py.
__version__ = '2.0.5'

__author__ = "Kyle S. Smith"

__doc__ = """\

API
======

Central class
-------------

.. autosummary::
   :toctree: .
   
   cfDNA


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