from __future__ import absolute_import, division, print_function
import os.path as op
import numpy as np
import pandas as pd
import numpy.testing as npt
import autoams as amps
from autoams.step4_insanemembrane import Insanity

data_path = op.join(amps.__path__[0], 'data')


def test_Insanity(pdb = 'all_molecules_combined.pdb'):
    '''Test INSANTIY class'''
    test = Insanity(pdb)
    test.MembraneSetup(l=['DOPE:15', 'DOPG:80', 'CDL1:5'], a=.5)
    test.PeriodicBoundarySetup(pbc='rectangular', x=22, y=22, z=22)
    test.SolventSetup(sol="PW", salt=0.15)
    test.ProteinSetup(dm="-9.9")
    assert test.command == "-l CDL1:5 -a 0.5 -pbc rectangular -x 22 -y 22 -z 22 -sol PW -salt 0.15 -dm -9.9"
