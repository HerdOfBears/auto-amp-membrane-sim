from __future__ import absolute_import, division, print_function
import os.path as op
import numpy as np
import pandas as pd
import numpy.testing as npt
import autoams as amps

data_path = op.join(amps.__path__[0], 'data')


def test_PeptideGrid(pdb = 'testpdb.pdb'):
    '''Test peptide grid object'''
    model = amps.PeptideGrid(pdb)    #create the object that will generate the grid of peptides
    model.GenerateGrid()        #generate the grid with default settings
    model.SaveGrid()            #save the 3d coordinates to a pdb file renamed to all_molecules_combined.pdb
    assert op.isfile('all_molecules_combined.pdb')
