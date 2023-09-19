from __future__ import absolute_import, division, print_function
import numpy as np
import pandas as pd
from .due import due, Doi

__all__ = [
    "simple_function"
]


# Use duecredit (duecredit.org) to provide a citation to relevant work to
# be cited. This does nothing, unless the user has duecredit installed,
# And calls this with duecredit (as in `python -m duecredit script.py`):
due.cite(Doi(""),
         description="Package for automated support of AMP and membrane simulations",
         tags=[""],
         path='autoams')


def simple_function(x, y):
    return x+y