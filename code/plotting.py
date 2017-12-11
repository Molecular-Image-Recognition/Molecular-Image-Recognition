#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 14:52:36 2017

@author: mattjohnson
"""
from lines_to_graph import lines_to_graph
from line import normalize,get_hough_lines
from skimage.feature import canny
from skimage.io import imread
from copy import deepcopy

class plotWrap:
    
    def __init__(self,path,params):
        self.path = path
        self.params = params
        
    def getPlottingValues(self):
        path = self.path
        params = self.params
        
        lines = get_hough_lines(canny(imread(path, as_grey=True)))
        lines_final = deepcopy(lines)
        imdim = max([x.length for x in lines])
        
        lines = normalize(lines,imdim)
        
        lines1 = lines_to_graph(lines, [0.2,0.4,1.0,0.2,0.5,0.5,0.5,0.2])
        
        imdim *= max([x.length for x in lines1])
        #imdim = np.median(np.array([x.length for x in lines]))
        
        # normalize all the lines
        lines_final = normalize(lines_final,imdim)
        #    min_dist_merge = params[0]
        #    min_angle_merge = params[1] 
        #    min_width_merge = params[2]
        #    split_tol = params[3]
        #    min_dist_bond = params[4]
        #    max_dist_bond = params[5]
        #    max_angle_bond = params[6]
        #    node_radius = params[7]
        lines = lines_to_graph(lines_final, params)
        
        xs = [[line.pts[0].x,line.pts[1].x] for line in lines]
        ys = [[line.pts[0].y,line.pts[1].y] for line in lines]
        
        return xs,ys