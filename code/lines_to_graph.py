
from line import *
from numba import jit

'''
INPUT
lines - a list of LineSegments and params is a list of the parameters
params - a list of the parameters which includes
    dist_merge - the minimum distance below which line segments are merged
    angle_merge - angle below which line segments are merged
'''
#@jit
def lines_to_graph(lines, params):
    
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
            order = max(line1.order,line2.order)
            dist, angle, width = line1.getDifference(line2)
            if dist < min_dist_merge and angle < min_angle_merge and width < min_width_merge:
                merged = combineLines([line1, line2])
                merged.order = order
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
            brokenlengths = []
            if len(broken) != 2:
                for a in broken:
                    brokenlengths.append(a.length)
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

    plotLines(lines)
    # deal with high order bonds
    i = 0
    while i < len(lines) - 1:
        j = i + 1
        while j < len(lines):
            
            line1 = lines[i]
            line2 = lines[j]
            
            dist, angle, width = line1.getDifference(line2)
            proj = line1.getProjOverlap(line2)
            
            if min_dist_bond < dist and dist < max_dist_bond and angle < max_angle_bond and proj > 0.2:
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
    
    iters = 0
    i = 0
    while i < len(lines) - 1:
        j = i + 1
        while j < len(lines):
            iters += 1
            if iters > 1000: #give up
                return lines
            line1 = lines[i]
            line2 = lines[j]
            
            a, b = line1.pts
            c, d = line2.pts
            L1c = line1.getShortestDistToPoint(c)
            L1d = line1.getShortestDistToPoint(d)
            L2a = line2.getShortestDistToPoint(a)
            L2b = line2.getShortestDistToPoint(b)

            if L1c > 1e-6 and L1c < node_radius:
                newlines = line1.pointSplit(c)
                lines[i] = newlines[0]
                lines.append(newlines[1])
                i = -1
                break
            elif L1d > 1e-6 and L1d < node_radius:
                newlines = line1.pointSplit(d)
                lines[i] = newlines[0]
                lines.append(newlines[1])
                i = -1
                break
            elif L2a > 1e-6 and L2a < node_radius:
                newlines = line2.pointSplit(a)
                lines[j] = newlines[0]
                lines.append(newlines[1])
            elif L2b > 1e-6 and L2b < node_radius:
                newlines = line2.pointSplit(b)
                lines[j] = newlines[0]
                lines.append(newlines[1])
            else:
                j += 1

        i += 1

    lines = combinePoints(lines, node_radius)
    
    for i,line1 in enumerate(lines):
        for j,line2 in enumerate(lines):
            if line1 == line2:
                continue
            elif line1.pts[0] == line2.pts[0] and line1.pts[1] == line2.pts[1]:
                lines[i] = line2
            elif line1.pts[1] == line2.pts[0] and line1.pts[0] == line2.pts[1]:
                lines[i] = line2
    
    lineList = []
    for line in lines:
        lineList.append((line,))
    
    lineList = list(set(lineList))
    
    lines = []
    for line in lineList:
        lines.append(line[0])
    
    plotLines(lines)
    
    return lines