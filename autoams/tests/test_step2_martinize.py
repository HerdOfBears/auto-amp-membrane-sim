from __future__ import absolute_import, division, print_function
import os.path as op
import numpy as np
import numpy.testing as npt
import autoams as amps
import logging

current_path = op.dirname(op.realpath(__file__))
data_path = op.join(amps.__path__[0], 'data')

def test_CoarseGrainer():
    
    cg = amps.step2_coarse_grain.CoarseGrainer()

    # load the test configuration file
    configs = cg.load_config_file(op.join(current_path, 'test_configuration_file.yml'))
    cg.set_configs(configs)
    assert cg.configs["-ff"] == "martini22p"

    # run the coarse graining
    cg.coarse_grain(cg.configs)