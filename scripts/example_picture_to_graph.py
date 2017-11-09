import os
from image_process.hough_transform import *
from regression.lines_to_graph import *
from molecule.adjacency import *

# image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'hexagon.JPG')

image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'Test_Set_1', 'PNGs',
                     'C(C)C(CCCC)(C)C.png')

lines = get_hough_lines(image)
imdim = max(get_image_size(image))

# normalize all the lines
for l in lines:

adj = lines_to_graph(lines, [])

atomnames = ['C']*adj.shape[0]
molecule = Adjacency(adj, atomnames)
molecule.addHydrogens()
rmgmol = molecule.toRMGmol()

