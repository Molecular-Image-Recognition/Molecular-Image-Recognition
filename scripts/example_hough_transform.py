from skimage import io
import os
from images.hough_transform import hough_transform_plots

image = io.imread(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'hexagon.JPG'),
                  as_grey=True)

hough_transform_plots(image, use_probabilistic=True)
