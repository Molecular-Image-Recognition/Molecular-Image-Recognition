import os
os.chdir('..')
from image_process.hough_transform import get_hough_lines
import matplotlib.pyplot as plt


#image_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'Test_Set_1', 'PNGs',
 #                         'C(C)C(CCCC)(C)C.png')
image_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test',
                          'hexagon.JPG')
line_segment_list = get_hough_lines(image_path)

l1 = line_segment_list[0]
l2 = line_segment_list[10]


print l1.pts
print l2.pts
print l1.getDifference(l2)
print l1.m
print l1.b

plt.figure()
plt.plot([l1.pts[0].x, l1.pts[1].x], [l1.pts[0].y, l1.pts[1].y])
plt.plot([l2.pts[0].x, l2.pts[1].x], [l2.pts[0].y, l2.pts[1].y])

deltas = []
thetas = []

for line_segs in line_segment_list:
    if line_segs == l1:
        continue

    diffs = l1.getDifference(line_segs)
    deltas.append(diffs[0])
    thetas.append(diffs[1])

plt.figure()
plt.plot(deltas, thetas, 'o')
plt.show()


