import os
os.chdir('..')
from image_process.hough_transform import hough_transform_plots, get_hough_lines

# image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'hexagon.JPG')

image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'Test_Set_1', 'PNGs',
                     'C(C)C(CCCC)(C)C.png')
image = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'hexagon.JPG')
hough_transform_plots(image)

print get_hough_lines(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test',
                                   'Test_Set_1', 'PNGs', 'C(C)C(CCCC)(C)C.png'))
