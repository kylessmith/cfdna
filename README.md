# Next Generation Sequencing fragment manipulation

[![Build Status](https://travis-ci.org/kylessmith/NGSfragments.svg?branch=master)](https://travis-ci.org/kylessmith/NGSfragments) [![PyPI version](https://badge.fury.io/py/NGSfragments.svg)](https://badge.fury.io/py/NGSfragments)

This is a Python package for easy and efficient manipulation
of Next Generation Sequencing reads.


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
    pip install NGSfragments
```

## Usage

```python
from NGSfragments import fragments
import numpy as np

# Create data
frags = fragments("test_bam.bam", verbose=True)
```
