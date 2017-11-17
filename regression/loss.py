#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 14:52:36 2017

@author: mattjohnson
"""
import os
os.chdir('..')

from rmgpy.molecule.molecule import Molecule
from molecule.adjacency import Adjacency
from lines.line import getAdjMatrix, LineSegment
from image_processing import get_hough_lines
from lines_to_graph import lines_to_graph

def getLoss(files,params):
    """
    calculates the loss for images in files and the set of parameters params
    """
    loss = 0.0
    for fname in files:
        if test(fname,params):
            pass
        else:
            loss += 1.0
    return loss

def test(fname,params):
    """
    checks if the given parameters allow proper identification of the image in fname
    assumes fname is of the form '/any/any/any/smilesString.any'
    for example /stuff/CC[O].thing 
    """
    bname = os.path.basename(os.path.normpath(fname)) #get rid of all but the final name
    smiles = bname.split('.')[0]
    m = Molecule().fromSMILES(smiles) #make actual molecule
    
    lines = get_hough_lines(fname)
    imdim = max([x.length for x in lines])

    # normalize all the lines
    for i,l in enumerate(lines):
        pt1 = l.pts[0]
        pt2 = l.pts[1]
        pt1.x /= imdim
        pt2.x /= imdim
        pt1.y /= imdim
        pt2.y /= imdim
        lines[i] = LineSegment([pt1,pt2])
    
    lines = lines_to_graph(lines, params)
    
    node_radius = params[6]
    
    adj = getAdjMatrix(lines,node_radius)
    atomnames = ['C']*adj.shape[0]
    molecule = Adjacency(adj, atomnames)
    molecule.addHydrogens()
    rmgmol = molecule.toRMGmol()
    
    return m.isIsomorphic(rmgmol)