from __future__ import absolute_import, division, print_function
import os.path as op
import numpy as np
import pandas as pd
import numpy.testing as npt
import autoams as amps

data_path = op.join(amps.__path__[0], 'data')


def test_simple_function_add():
    assert amps.simple_function_add(1, 2) == 3
    assert amps.simple_function_add("a", "b") == "ab"

def test_simple_function_new_branch():
    assert amps.simple_function_new_branch(1, 2) == 1
    assert amps.simple_function_new_branch(2,2)  == 4