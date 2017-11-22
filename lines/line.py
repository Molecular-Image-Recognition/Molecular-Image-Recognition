
import numpy as np

class Point:
    def __init__(self,x,y):
        self.x = float(x)
        self.y = float(y)
        
    def __repr__(self):
        return '({0},{1})'.format(self.x,self.y)
    
    def getDistance(self,pt):
        return np.sqrt((self.x-pt.x)**2+(self.y-pt.y)**2)

class LineSegment:
    def __init__(self,pts,order=1):
        self.pts = pts
        if abs(pts[1].x-pts[0].x) > 1e-15:
            self.m = (pts[1].y-pts[0].y)/(pts[1].x-pts[0].x)
        else:
            self.m = np.inf
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
        boo1 = pt.y-self.m*pt.x-self.b < 1e-4
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
        return min(dists)
    
    def getArea(self,L):
        pts = self.pts+L.pts
        return getArea(pts)
    

    def getDifference(self,L):
        """
        find important difference information between two lines
        """
        A = self.getArea(L)
        length = combineLines([self,L]).length
        width = A/length
        thta = abs(self.theta-L.theta)
        thta = min(thta,np.pi-thta)
        return [self.getMinimumDist(L),thta,width]
    
    def getIntersection(self,L):
        """
        find the point at which two lines intersect
        """
        if self.m == L.m:
            return None
        else:
            x = (L.b-self.b)/(self.m-L.m)
            y = self.m*x+self.b
            pt = Point(x,y)
            if self.pointIn(pt) and L.pointIn(pt):
                return Point(x,y)
            else:
                return None
    
    def breakAtIntersection(self,L):
        """
        break the two lines at their itersection into four lines
        """
        pt = self.getIntersection(L)
        if pt is None:
            return [self,L]
        notok = (self.pts[0].x == pt.x and self.pts[0].y == pt.y) or (self.pts[1].x == pt.x and self.pts[1].y == pt.y)
        if not notok and self.pointIn(pt):
            lines = [LineSegment([self.pts[0],pt]),LineSegment([self.pts[1],pt]),LineSegment([L.pts[0],pt]),LineSegment([L.pts[1],pt])]
            return lines
        else:
            return [self,L]
            
    
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
#    if m >= 0:
#        xmaxy = m*xmax+b
#        if xmaxy > ymax:
#            pt1 = Point(xmax,xmaxy)
#        else:
#            pt1 = Point((ymax-b)/m,ymax)
#        xminy = m*xmin+b
#        if xminy < ymin:
#            pt2 = Point(xmin,xminy)
#        else:
#            pt2 = Point((ymin-b)/m,ymin)
#    else:
#        xminy = m*xmin+b
#        if xminy > ymax:
#            pt1 = Point(xmin,xminy)
#        else:
#            pt1 = Point((ymax-b)/m,ymax)
#        xmaxy = m*xmax+b
#        if xmaxy < ymin:
#            pt2 = Point(xmax,xmaxy)
#        else:
#            pt2 = Point((ymin-b)/m,ymin)
    
#    f = lambda x: m*x + b
#    finv = lambda y: (y-b)/m
#    
#    if ymin < f(xmin) < ymax:
#        pt1 = Point(xmin, f(xmin))
#    else:
#        yp = ymin if m > 0 else ymax
#        pt1 = Point(finv(yp), yp)
#        
#    if ymin < f(xmax) < ymax:
#        pt2 = Point(xmax, f(xmax))
#    else:
#        yp = ymax if m > 0 else ymin
#        pt2 = Point(finv(yp), yp)
    s = 0.0
    pts = []
    for line in lines:
        pts += line.pts
    for pt1 in pts:
        for pt2 in pts:
            val = pt1.getDistance(pt2)
            if val > s:
                outpts =[pt1,pt2]
                s = val
    return LineSegment(outpts)

def combinePoints(lines,atol):
    """
    adjusts the points so lines close to itersecting intersect
    """
    pts = []
    for line in lines:
        pts += line.pts
    
    identicals = []
    for i,pt1 in enumerate(pts):
        for j,pt2 in enumerate(pts):
            if i != j:
                dist = pt1.getDistance(pt2)
                if dist < atol and pt1 != pt2:
                    identicals.append({i,j})
                    
    for i,s in enumerate(identicals):
        for j,s2 in enumerate(identicals):
            if s & s2 != set():
                identicals[i] = identicals[i] | s2
                identicals[j] = identicals[j] | s
    
    identicals = [tuple(item) for item in identicals]
    identicals = list(set(identicals))
    
    for s in identicals:
        pt = pts[s[0]]
        for i in xrange(1,len(s)):
            pts[s[i]] = pt
    
    lines = [LineSegment([pt,pts[i+1]]) for i,pt in enumerate(pts) if i%2 == 0]
    delinds = [] #delete lines for which points are identical
    for i,line in enumerate(lines):
        if line.pts[0] == line.pts[1]:
            delinds.append(i)
    
    for i in reversed(delinds):
        del lines[i]
        
    return lines

def getAdjMatrix(lines,atol):
    """
    converts a set of lines with point objects that have been made identical if close
    into an adjacency matrix, atol is the point closeness tolerance
    """
    
    lines = __combinePoints(lines,atol)
    
    pts = []
    for line in lines:
        pts += lines.pts
    
    adjMat = np.zeros((len(pts),len(pts)))
    
    for i,pt in enumerate(pts):
        for line in lines:
            if line.pts[0] == pt:
                adjMat[i,pts.index(line.pts[0])] = line.order
                adjMat[pts.index(line.pts[0]),i] = line.order
            elif line.pts[1] == pt:
                adjMat[i,pts.index(line.pts[1])] = line.order
                adjMat[pts.index(line.pts[1]),i] = line.order
    
    return adjMat

def getArea(pts):
    """
    area calculated from verticies of 2-D shape based on the shoelace formula
    """
    N = len(pts)
    s1 = pts[N-1].x*pts[0].y
    s2 = pts[0].x*pts[N-1].y
    s1 += sum([pts[i].x*pts[i+1].y for i in xrange(0,N-1)])
    s2 += sum([pts[i+1].x*pts[i].y for i in xrange(0,N-1)])
    
    return .5*abs(s1-s2)

    
