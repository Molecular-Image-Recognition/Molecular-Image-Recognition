#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 18:14:27 2017

@author: mattjohnson
"""
import os
from line import get_hough_lines
from lines_to_graph import *
from skimage.feature import canny
from skimage.io import imread
from loss import test
from copy import deepcopy
from skimage.filters import threshold_otsu
import numpy as np

class Dataset:
    
    def __init__(self,directory):
        paths = os.listdir(directory)
        paths = [os.path.join(directory,path) for path in paths]
        linesList = []
        smiles = []
        for path in paths:
            bname = os.path.basename(os.path.normpath(path)) #get rid of all but the final name
            smile = bname.split('.')[0]
            smiles.append(smile)
            
            image = imread(path, as_grey=True)
            thresh = threshold_otsu(image)
            image = image > thresh
            image = np.invert(image)
            
            linesList.append(get_hough_lines(image))
        
        print smiles 
        self.smiles = smiles

#        for j,lines in enumerate(linesList):
#        
#            imdim = max([x.length for x in lines])
#            #imdim = np.median(np.array([x.length for x in lines]))
#            
#            lines_copy = deepcopy(lines)
#            # normalize all the lines
#            for i,l in enumerate(lines):
#                pt1 = l.pts[0]
#                pt2 = l.pts[1]
#                pt1.x /= imdim
#                pt2.x /= imdim
#                pt1.y /= imdim
#                pt2.y /= imdim
#                lines[i] = LineSegment([pt1,pt2])
#                
#                
#            lines1 = lines_to_graph(lines, [0.2,0.4,1.0,0.2,0.5,0.5,0.5,0.2])
#            
#            imdim *= max([x.length for x in lines1])
#            
#            for i,l in enumerate(lines_copy):
#                pt1 = l.pts[0]
#                pt2 = l.pts[1]
#                pt1.x /= imdim
#                pt2.x /= imdim
#                pt1.y /= imdim
#                pt2.y /= imdim
#                lines_copy[i] = LineSegment([pt1,pt2])
#            
#            linesList[j] = lines_copy
        
        
        self.linesList = linesList
        
        for lines in linesList:
            print len(lines)
        
    
    def getLoss(self,params):
        """
        calculates the loss for images in dataset and the set of parameters params
        """
        loss = 0.0
        linesList = self.linesList
        smiles = self.smiles
        for i,lines in enumerate(linesList):
            if test(lines,smiles[i],params):
                pass
            else:
                loss += 1.0
        
        return loss
    
