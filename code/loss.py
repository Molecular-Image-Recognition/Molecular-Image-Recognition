#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 14:52:36 2017

@author: mattjohnson
"""
import os
import numpy as np
from rmgpy.molecule.molecule import Molecule
from adjacency import Adjacency
from line import *
from line import get_hough_lines
from lines_to_graph import lines_to_graph
from skimage.feature import canny
from skimage.io import imread
from copy import deepcopy
from numba import jit

def test(lines,smiles,params):
    """
    checks if the given parameters allow proper identification of the image
    based on comparison with the smiles and the use of the set of parameters
    """

    m = Molecule().fromSMILES(smiles) #make actual molecule

    lines = lines_to_graph(lines, params)
    
    
    adj = getAdjMatrix(lines)
    atomnames = ['C']*adj.shape[0]
    try:
        molecule = Adjacency(adj, atomnames)
    
        molecule.addHydrogens()
        rmgmol = molecule.toRMGmol()
    except:
        print 'failed to generate molecule'
        print smiles
        plt.figure()
        plotLines(lines)
        return False
    boo = m.isIsomorphic(rmgmol)
    if not boo:
        print smiles
        plt.figure()
        plotLines(lines)
    return boo