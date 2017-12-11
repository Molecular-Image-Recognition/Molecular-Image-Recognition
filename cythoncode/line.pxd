#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 20:19:18 2017

@author: mattjohnson
"""
cimport numpy

cdef class Point(object):
    cdef public float x
    cdef public float y
    cpdef getDistance(self,object pt)
    
cdef class LineSegment(object):
    cdef public list pts
    cdef public float m
    cdef public float b
    cdef public float theta
    cdef public float xmin
    cdef public float xmax
    cdef public float length
    cdef public int order
    cpdef pointIn(self,object pt)
    cpdef getArea(self,object L)
    cpdef getShortestDistToPoint(self,object pt)
    cpdef getMinimumDist(self,object L)
    cpdef getDifference(self, object L)
    cpdef getIntersection(self,object L)
    cpdef breakAtIntersection(self,object L)
    cpdef extendToItersection(self,object L)
    
cpdef combineLines(list lines)

cpdef combinePoints(list lines, float atol)

cpdef getAdjMatrix(list lines)

cpdef getArea(list pts)

cpdef get_hough_lines(numpy.ndarray edgesn)