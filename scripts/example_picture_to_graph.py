import os
os.chdir('..')
from image_process.hough_transform import *
from regression.lines_to_graph import *
from molecule.adjacency import *
from lines.line import *
import matplotlib.pyplot as plt
# image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'hexagon.JPG')

#image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'Test_Set_1', 'PNGs',
                   #  'C(C)C(CCCC)(C)C.png')

image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test','square.jpg' )

lines = get_hough_lines(image)
#imsize = min(get_image_size(image))

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
    
    
#    min_dist_merge = params[0]
#    min_angle_merge = params[1]
#    split_tol = params[2]
#    min_dist_bond = params[3]
#    max_dist_bond = params[4]
#    max_angle_bond = params[5]
#    node_radius = params[6]
lines = lines_to_graph(lines, [1.0,0.5,0.1,0.5,0.5,0.5,.1])

for line in lines:
    plt.plot((line.pts[0].x, line.pts[1].x), (line.pts[0].y, line.pts[1].y))
plt.xlim((0, get_image_size(image)[1]/imdim))
plt.ylim(( get_image_size(image)[0]/imdim,0))

plt.title('Probabilistic Hough')

adj = getAdjMatrix(lines,node_radius)
atomnames = ['C']*adj.shape[0]
molecule = Adjacency(adj, atomnames)
molecule.addHydrogens()
rmgmol = molecule.toRMGmol()

