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
#    min_dist_merge = params[0]
#    min_angle_merge = params[1] 
#    min_width_merge = params[2]
#    split_tol = params[3]
#    min_dist_bond = params[4]
#    max_dist_bond = params[5]
#    max_angle_bond = params[6]
#    node_radius = params[7]

directory = '/Users/mattjohnson/RMGCODE/Molecular-Image-Recognition/data/hand_drawn'

#testdirectory = '/Users/mattjohnson/RMGCODE/Molecular-Image-Recognition/data/hand_drawn_test'

#directory = '/Users/mattjohnson/RMGCODE/Molecular-Image-Recognition/data/machine_drawn_all'

#testdirectory = '/Users/mattjohnson/RMGCODE/Molecular-Image-Recognition/data/hand_drawn_test'
#bounds = [(0,1),(0,np.pi/2.0),(0,1),(0,0.5),(0,1),(0,1),(0,np.pi/2.0),(0,1)]

bounds = [(0,.16),(.3,.5),(.1,.3),(0,.4),(0,.1),(.3,.4),(.2,.4),(.15,.25)]


vals = []
for i in xrange(1):
    dset = Dataset(directory)

    print 'loaded data'
#print len(dset.linesList)
#dset.linesList = dset.linesList[0:5]
#dset.smiles = dset.smiles[0:5]

#start = time.time()
    
    val = dset.getLoss([ 0.1393534 ,  0.38026248,  0.2576827 ,  0.16425988,  0.07834396,
        0.32466892,  0.26150276,  0.21460726])
#    print val
#    vals.append(val)
#
#print min(vals)
#print max(vals)
#print sum(vals)/10.0
#out = differential_evolution(dset.getLoss,bounds,popsize=3,maxiter=10,tol=1e-16)

#print out

#out = differential_evolution(dset.getLoss,bounds,popsize=5,maxiters=8)