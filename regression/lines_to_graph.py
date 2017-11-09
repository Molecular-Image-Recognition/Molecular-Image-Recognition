import numpy as np
from image_process.hough_transform import *
from lines.line import *

'''
INPUT
lines - a list of LineSegments and params is a list of the parameters
params - a list of the parameters which includes
    dist_merge - the minimum distance below which line segments are merged
    angle_merge - angle below which line segments are merged
'''
def lines_to_graph(lines, params):
    # define paramters
    min_dist_merge = params[0]
    min_angle_merge = params[1]
    split_tol = params[2]
    min_dist_bond = params[3]
    max_dist_bond = params[4]
    max_angle_bond = params[5]
    node_radius = params[6]
    
    i=0
	# merge lines
    while i < len(lines) - 1:
        j = i + 1
        while j < len(lines):
            print(i)
            print(j)
            didmerge = False

            line1 = lines[i]
            line2 = lines[j]
            
            dist, angle = line1.getDifference(line2)
            if dist < min_dist_merge and angle < min_angle_merge:
                merged = combineLines([line1, line2])
                lines[i] = merged
                del lines[j]
                didmerge = True

            if not didmerge:
                j += 1
        i = i + 1
    return lines
    i=0
    # deal with intersections
    while i < len(lines) - 1:
        j = i + 1
        while j < len(lines):
            print(i)
            print(j)
            didbreak = False
            line1 = lines[i]
            line2 = lines[j]

            broken = line1.breakAtIntersection(line2)

            if len(broken) != 2:
                didbreak = True
                brokenlengths = [a.length for a in broken]
                brokenlengths[0] /= line1.length
                brokenlengths[1] /= line1.length
                brokenlengths[2] /= line2.length
                brokenlengths[3] /= line2.length

                for k in [3,2,1,0]:
                    if brokenlengths[k] < split_tol:
                        del broken[k]

            if didbreak:
                lines[i] = broken[0]
                lines[j] = broken[1:]

                j += len(broken) - 1
            else:      
                j += 1

	# join bonds together and return the adjacency matrix generated
        
    return lines