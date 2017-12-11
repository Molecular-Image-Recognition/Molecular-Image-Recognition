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
from line import getAdjMatrix, LineSegment
from hough_transform import get_hough_lines
from lines_to_graph import lines_to_graph
from skimage.feature import canny
from skimage.io import imread
from copy import deepcopy
import cython

cpdef test(list lines,str smiles,list params):
    """
    checks if the given parameters allow proper identification of the image
    based on comparison with the smiles and the use of the set of parameters
    """
    cdef object m,pt1,pt2,x,adj,molecule,rmgmol
    cdef float imdim
    cdef list lines_final,lines1,atomnames

    m = Molecule().fromSMILES(smiles) #make actual molecule
    
    imdim = max([x.length for x in lines])
    
    lines_final = deepcopy(lines)
    
    imdim = max([x.length for x in lines])
    #imdim = np.median(np.array([x.length for x in lines]))
    
    # normalize all the lines
    for i,l in enumerate(lines):
        pt1 = l.pts[0]
        pt2 = l.pts[1]
        pt1.x /= imdim
        pt2.x /= imdim
        pt1.y /= imdim
        pt2.y /= imdim
        lines[i] = LineSegment([pt1,pt2])
        
        
    lines1 = lines_to_graph(lines, [0.2,0.4,1.0,0.2,0.5,0.5,0.5,0.2])
    
    imdim *= max([x.length for x in lines1])
    
    # normalize all the lines
    for i,l in enumerate(lines_final):
        pt1 = l.pts[0]
        pt2 = l.pts[1]
        pt1.x /= imdim
        pt2.x /= imdim
        pt1.y /= imdim
        pt2.y /= imdim
        lines_final[i] = LineSegment([pt1,pt2])
    
    lines = lines_to_graph(lines_final, params)

    adj = getAdjMatrix(lines)
    atomnames = ['C']*adj.shape[0]
    
    try:
        molecule = Adjacency(adj, atomnames)
        molecule.addHydrogens()
        rmgmol = molecule.toRMGmol()
    except:
        return False
    
    return m.isIsomorphic(rmgmol)