import numpy as np
from skimage.transform import hough_line, hough_line_peaks, probabilistic_hough_line
from skimage.feature import canny
from skimage.io import imread
from lines.line import Point, LineSegment
import matplotlib.pyplot as plt


# Code below modified from code at  http://scikit-image.org/docs/dev/auto_examples/edges/plot_line_hough_transform.html

def hough_transform_plots(image_path, use_probabilistic=True, use_canny=True):

    edges = _load_image(image_path, use_canny)

    image = imread(image_path)

    h, theta, d = hough_line(edges)
    lines = probabilistic_hough_line(edges)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), )
    ax = axes.ravel()

    ax[0].imshow(image, cmap=plt.gray())
    ax[0].set_title('Input image')

    ax[1].imshow(edges, cmap=plt.gray())
    ax[1].set_title('Canny edges')

    if use_probabilistic:
        ax[2].imshow(edges * 0)
        for line in lines:
            p0, p1 = line
            ax[2].plot((p0[0], p1[0]), (p0[1], p1[1]))
        ax[2].set_xlim((0, image.shape[1]))
        ax[2].set_ylim((image.shape[0], 0))
        ax[2].set_title('Probabilistic Hough')

    else:
        ax[2].imshow(image)
        for _, angle, dist in zip(*hough_line_peaks(h, theta, d)):
            y0 = (dist - 0 * np.cos(angle)) / np.sin(angle)
            y1 = (dist - image.shape[1] * np.cos(angle)) / np.sin(angle)
            ax[2].plot((0, image.shape[1]), (y0, y1), '-r')
        ax[2].set_xlim((0, image.shape[1]))
        ax[2].set_ylim((image.shape[0], 0))
        ax[2].set_axis_off()
        ax[2].set_title('Detected lines')

    for a in ax:
        a.set_axis_off()
        a.set_adjustable('box-forced')

    plt.tight_layout()
    plt.show()


def get_hough_lines(image_path, use_canny=True):
    edges = _load_image(image_path, use_canny)
    lines = probabilistic_hough_line(edges)

    line_segment_list = []

    for line_seg in lines:
        point_1 = Point(*line_seg[0])
        point_2 = Point(*line_seg[1])
        new_line = LineSegment([point_1, point_2])
        line_segment_list.append(new_line)

    return line_segment_list


def _load_image(image_path, use_canny):
    image = imread(image_path, as_grey=True)

    if use_canny:
        edges = canny(image)
    else:
        edges = image
    return edges


def get_image_size(image_path):
    image = imread(image_path, as_grey=True)
    return image.shape
