from image_process.hough_transform import get_hough_lines
from lines.line import Point, LineSegment
import matplotlib.pyplot as plt
import os

image_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data', 'hough_test', 'Test_Set_1', 'PNGs',
                          'C(C)C(CCCC)(C)C.png')

lines = get_hough_lines(image_path)

line_segment_list = []

for line_seg in lines:
    point_1 = Point(*line_seg[0])
    point_2 = Point(*line_seg[1])
    new_line = LineSegment([point_1, point_2])
    line_segment_list.append(new_line)

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
