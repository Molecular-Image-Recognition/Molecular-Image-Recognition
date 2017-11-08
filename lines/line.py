
import numpy as np

class Point:
    def __init__(self,x,y):
        self.x = x
        self.y = y
        
    def __repr__(self):
        return '({0},{1})'.format(self.x,self.y)
    
    def getDistance(self,pt):
        return np.sqrt((self.x-pt.x)**2+(self.y-pt.y)**2)

class LineSegment:
    def __init__(self,pts,order=1):
        self.pts = pts
        self.m = (pts[1].y-pts[0].y)/(pts[1].x-pts[0].x)
        self.b = self.pts[0].y-self.m*self.pts[0].x
        self.theta = np.arctan(self.m)
        xs = [pt.x for pt in pts]
        #ys = [pt.y for pt in pts]
        self.xmin = min(xs)
        self.xmax = max(xs)
        #self.ymin = min(ys)
        #self.ymax = max(ys)
        self.length = pts[0].getDistance(pts[1])
        self.order = order #bond order
        
    def pointIn(self,pt):
        """
        is pt within the line segement
        """
        boo1 = pt.y-self.m*pt.x-self.b < 1e-8
        boo2 = pt.x >= self.xmin and pt.x <= self.xmax
        return boo1 and boo2
    
    def getShortestDistToPoint(self,pt):
        """
        The shortest distance between the point pt and the overall line (not segemented)
        """
        mp = -1.0/self.m
        bp = pt.y-mp*pt.x
        xstar = (bp-self.b)/(self.m-mp)
        ystar = self.m*xstar+self.b
        ptstar = Point(xstar,ystar)
        boo =  self.pointIn(ptstar)
        dist = ptstar.getDistance(pt)
        if boo:
            return dist
        else:
            return np.inf
    
    def getMinimumDist(self,L):
        """
        find the minimum distance between two line segements
        """
        dists = [self.getShortestDistToPoint(pt) for pt in L.pts]
        dists += [L.getShortestDistToPoint(pt) for pt in self.pts]
        dists += [L.pts[0].getDistance(pt) for pt in self.pts]
        dists += [L.pts[1].getDistance(pt) for pt in self.pts]
        print dists
        return min(dists)
    
    def getDifference(self,L):
        """
        find important difference information between two lines
        """
        return [self.getMinimumDist(L),abs(self.theta-L.theta)]
    
    def getIntersection(self,L):
        """
        find the point at which two lines intersect
        """
        if self.m == L.m:
            return None
        else:
            x = (L.b-self.b)/(self.m-L.m)
            y = self.m*x+self.b
            return Point(x,y)
    
    def breakAtItersection(self,L):
        """
        break the two lines at their itersection into four lines
        """
        pt = self.getIntersection(L)
        if pt and self.pointIn(pt):
            lines = [LineSegment(self.pts[0],pt),LineSegment(self.pts[1],pt),LineSegment(L.pts[1],pt),LineSegment(L.pts[1],pt)]
            return lines
        else:
            return [self]
    
    def extendToItersection(self,L):
        """
        extend two lines to the point at which they intersect
        """
        pt = self.getItersection(L)
        if pt:
            dist = [pt.getDistance(lpt) for lpt in self.pts]
            indself = dist.index(min(dist))
            dist = [pt.getDistance(lpt) for lpt in L.pts]
            indL = dist.index(min(dist))
            return [LineSegment(pt,self.pts[indself]),LineSegment(pt,L.pts[indL])]
        else:
            return [self,L]
        
def combineLines(lines):
    """
    combine the set of LineSegments in lines into one line segment
    """
    xs = []
    ys = []
    for line in lines:
        xs.append(line.pts[0].x)
        xs.append(line.pts[1].x)
        ys.append(line.pts[0].y)
        ys.append(line.pts[1].y)
    A = np.ones((len(ys),2))
    A[:,0] = np.array(xs)
    y = np.array(ys)
    m,b = np.linalg.lstsq(A,y)[0]
    xmax = max(xs)
    xmin = min(xs)
    ymax = max(ys)
    ymin = min(ys)
    if m >= 0:
        xmaxy = m*xmax+b
        if xmaxy > ymax:
            pt1 = Point(xmax,xmaxy)
        else:
            pt1 = Point((ymax-b)/m,ymax)
        xminy = m*xmin+b
        if xminy < ymin:
            pt2 = Point(xmin,xminy)
        else:
            pt2 = Point((ymin-b)/m,ymin)
    else:
        xminy = m*xmin+b
        if xminy > ymax:
            pt1 = Point(xmin,xminy)
        else:
            pt1 = Point((ymax-b)/m,ymax)
        xmaxy = m*xmax+b
        if xmaxy < ymin:
            pt2 = Point(xmax,xmaxy)
        else:
            pt2 = Point((ymin-b)/m,ymin)
    
    return Line([pt1,pt2])

def combinePoints(lines,atol):
    """
    adjusts the points so lines close to itersecting intersect
    """
    pts = []
    lineList = []
    for line in lines:
        pts += lines.pts
        lineList.append(line)
    
    for i,pt1 in enumerate(pts):
        for j,pt2 in enumerate(pts):
            if i != j:
                dist = pt1.getDistance(pt2)
                if dist < atol:
                    ind = i // 2
                    ptind = i % 2
                    line = lineList[ind]
                    if ptind == 0:
                        newLine = LineSegment([pt2,line.pts[1]])
                    else:
                        newLine = LineSegment([pt2,line.pts[0]])
                    lineList[i] = newLine
    
    return lineList
