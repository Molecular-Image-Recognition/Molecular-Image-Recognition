#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 18:14:27 2017

@author: mattjohnson
"""
import os
from line import get_hough_lines
from skimage.feature import canny
from skimage.io import imread
from loss import test
import cython

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
            
            linesList.append(get_hough_lines(canny(imread(path, as_grey=True))))
            
        self.linesList = linesList
        self.smiles = smiles
    
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
    
