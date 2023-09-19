from __future__ import absolute_import, division, print_function
import os.path as op
import numpy as np
import pandas as pd
import numpy.testing as npt
import autoams as amps

data_path = op.join(amps.__path__[0], 'data')


def test_simple_function():
    assert amps.simple_function(1, 2) == 3
    assert amps.simple_function("a", "b") == "ab"