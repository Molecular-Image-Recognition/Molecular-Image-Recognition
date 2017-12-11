#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import os
import numpy as np
from lines_to_graph import *
from adjacency import *
from line import *
import matplotlib.pyplot as plt
from copy import deepcopy
from skimage.feature import canny
from skimage.io import imread

#image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'hexagon.JPG')
#image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'PDD.png')
#image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'Test_Set_1', 'PNGs',
#                     'C(C)C(CCCC)(C)C.png')

#image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'Test_Set_1', 'PNGs',
#                     'C(C)(CC)(CC)CCCCC.png')
image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test','square.jpg' )
#image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'Test_Set_1','double_bonds', 'C((CCC)CC)C.png')
#image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'Test_Set_1','double_bonds', 'C(C)CC.png')


lines = get_hough_lines(canny(imread(image, as_grey=True)))

lines_final = deepcopy(lines)
#for line in lines_final:
#    plt.plot((line.pts[0].x, line.pts[1].x), (line.pts[0].y, line.pts[1].y))
#plt.xlim((0, get_image_size(image)[1]))
#plt.ylim(( get_image_size(image)[0],0))
#plt.figure()
#imsize = min(get_image_size(image))

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
#imdim = np.median(np.array([x.length for x in lines]))

# normalize all the lines
for i,l in enumerate(lines_final):
    pt1 = l.pts[0]
    pt2 = l.pts[1]
    pt1.x /= imdim
    pt2.x /= imdim
    pt1.y /= imdim
    pt2.y /= imdim
    lines_final[i] = LineSegment([pt1,pt2])
    
#    min_dist_merge = params[0]
#    min_angle_merge = params[1] 
#    min_width_merge = params[2]
#    split_tol = params[3]
#    min_dist_bond = params[4]
#    max_dist_bond = params[5]
#    max_angle_bond = params[6]
#    node_radius = params[7]
lines = lines_to_graph(lines_final, [0.05,.4,.1,0.2,0.06,0.8,0.1,.2])
#lines = lines_to_graph(lines_final, [0.05,.4,.1,0.2,0.06,0.8,0.1,.2])
print [line.order for line in lines]
for line in lines:
    plt.plot((line.pts[0].x, line.pts[1].x), (line.pts[0].y, line.pts[1].y))
plt.xlim((0, get_image_size(image)[1]/imdim))
plt.ylim(( get_image_size(image)[0]/imdim,0))
plt.figure()

for line in lines:
    plt.plot((line.pts[0].x, line.pts[1].x), (line.pts[0].y, line.pts[1].y))
plt.xlim((0, get_image_size(image)[1]/imdim))
plt.ylim(( get_image_size(image)[0]/imdim,0))

plt.title('Probabilistic Hough')

plt.save()
adj = getAdjMatrix(lines)

try:
    molecule = Adjacency(adj, atomnames)
    molecule.addHydrogens()
    rmgmol = molecule.toRMGmol()
except:
     print 'failed to generate molecule'
    
print m.isIsomorphic(rmgmol)

