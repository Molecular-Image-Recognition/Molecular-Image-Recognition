#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 15:11:13 2017

@author: mattjohnson
"""
import numpy as np
from scipy.optimize import differential_evolution
from dataset import Dataset
import time
import cython
#    min_dist_merge = params[0]
#    min_angle_merge = params[1] 
#    min_width_merge = params[2]
#    split_tol = params[3]
#    min_dist_bond = params[4]
#    max_dist_bond = params[5]
#    max_angle_bond = params[6]
#    node_radius = params[7]

directory = '/Users/mattjohnson/RMGCODE/Molecular-Image-Recognition/data/hough_test/Test_Set_1/PNGs'

bounds = [(0,1),(0,np.pi/2.0),(0,1),(0,0.5),(0,1),(0,1),(0,np.pi/2.0),(0,1)]

dset = Dataset(directory)

#print 'loaded data'
#print len(dset.linesList)
dset.linesList = dset.linesList[0:4]
dset.smiles = dset.smiles[0:4]

start = time.time()

for i in range(1):

    dset.getLoss([b[1] for b in bounds])
    
end = time.time()

print end - start

#out = differential_evolution(dset.getLoss,bounds)