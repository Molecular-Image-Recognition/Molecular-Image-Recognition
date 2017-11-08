import cairosvg
import os

# MUST RUN IN PYTHON 3 and pip install cairosvg

file_dir = '../data/hough_test/Test_Set_1/'

svgs = os.listdir(os.path.join(file_dir, 'SVGs'))

for svg in svgs:
    name = svg.split('.svg')[0]
    cairosvg.svg2png(url=os.path.join(file_dir, 'SVGs', svg),
                     write_to=os.path.join(file_dir, 'PNGs', '{0}.png'.format(name)), dpi=600)
    # cairosvg.svg2pdf(url=os.path.join(file_dir, 'SVGs', svg),
                     # write_to=os.path.join(file_dir, 'PDFs', '{0}.pdf'.format(name)), dpi=600)