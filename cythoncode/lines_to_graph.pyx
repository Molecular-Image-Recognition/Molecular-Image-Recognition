import numpy as np
from line import *
import cython
cimport numpy

'''
INPUT
lines - a list of LineSegments and params is a list of the parameters
params - a list of the parameters which includes
    dist_merge - the minimum distance below which line segments are merged
    angle_merge - angle below which line segments are merged
'''

cpdef lines_to_graph(list lines, list params):
    
    cdef float min_dist_merge,min_angle_merge,min_width_merge,split_tol,min_dist_bond
    cdef float max_dist_bond,max_angle_bond,node_radius
    cdef int i,j,k
    cdef object line1,line2,merged
    cdef float dist,angle,width
    cdef list broken,brokenlengths
    
    # define paramters
    min_dist_merge = params[0]
    min_angle_merge = params[1]
    min_width_merge = params[2]
    split_tol = params[3]
    min_dist_bond = params[4]
    max_dist_bond = params[5]
    max_angle_bond = params[6]
    node_radius = params[7]

    i=0
	# merge lines
    while i < len(lines):
        j = i + 1
        while j < len(lines):
            didmerge = False

            line1 = lines[i]
            line2 = lines[j]
            
            dist, angle, width = line1.getDifference(line2)
            if dist < min_dist_merge and angle < min_angle_merge and width < min_width_merge:
                merged = combineLines([line1, line2])
                lines[i] = merged
                del lines[j]
                didmerge = True

            if didmerge:
                i = -1
                break
            else:
                j += 1
        i += 1
        
    # deal with intersections
    i = 0
    while i < len(lines) - 1:
        j = i + 1
        while j < len(lines):
            line1 = lines[i]
            line2 = lines[j]
            broken = line1.breakAtIntersection(line2)
            
            if len(broken) != 2:
                brokenlengths = [a.length for a in broken]
                brokenlengths[0] /= line1.length
                brokenlengths[1] /= line1.length
                brokenlengths[2] /= line2.length
                brokenlengths[3] /= line2.length

                for k in [3,2,1,0]:
                    if brokenlengths[k] < split_tol:
                        del broken[k]
                
                lines[i] = broken[0]
                lines[j] = broken[1]
                for k in reversed(range(2,len(broken))):
                    lines.insert(j+1, broken[k])

                j += len(broken) - 1
            else:      
                j += 1
        
        i += 1
        
    # deal with high order bonds
    i = 0
    while i < len(lines) - 1:
        j = i + 1
        while j < len(lines):
            
            line1 = lines[i]
            line2 = lines[j]
            
            dist, angle, width = line1.getDifference(line2)
            
            if min_dist_bond < dist and dist < max_dist_bond and angle < max_angle_bond:
                # remove the shorter of the lines and increment order
                if line1.length < line2.length:
                    lines[j].order += 1
                    del lines[i]
                else:
                    lines[i].order += 1
                    del lines[j]
                    
                i = -1
                break
            else:
                j += 1
        i += 1
                    

	# join bonds together and return the adjacency matrix generated
    lines = combinePoints(lines, node_radius)
        
    return lines