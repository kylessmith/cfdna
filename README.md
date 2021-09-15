# Cell-free DNA analysis toolkit

[![Build Status](https://travis-ci.org/kylessmith/cfdna.svg?branch=master)](https://travis-ci.org/kylessmith/cfdna) [![PyPI version](https://badge.fury.io/py/cfdna.svg)](https://badge.fury.io/py/cfdna)

This is a Python package for easy and efficient cell-free
DNA analysis.


## Install

If you dont already have numpy and scipy installed, it is best to download
`Anaconda`, a python distribution that has them included.  
```
    https://continuum.io/downloads
```

Dependencies can be installed by:

```
    pip install -r requirements.txt
```

PyPI install, presuming you have all its requirements installed:
```
    pip install cfdna
```

## Usage

```python
from cfdna import cfDNA
import numpy as np

# Create data
cf = cfDNA("test_bam.bam", verbose=True)
```
