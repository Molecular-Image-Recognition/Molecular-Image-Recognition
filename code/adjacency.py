import os
import numpy as np
from rmgpy.molecule.molecule import Molecule



class Adjacency:
    """
    class for holding molecular graph information and conversion to RMG objects
    """
    atomDict = {'H':1,'O':2,'S':2,'N':3,'C':4}
    
    def __init__(self, adjMatrix,
                 atomNames, unpaireds=None, lonepairs=None, charges=None):
        self.adjMatrix = adjMatrix
        self.atomNames = atomNames
        self.unpaireds = unpaireds
        self.lonepairs = lonepairs
        self.charges = charges
        self.AdjList = None
        self.RMGspc = None
        
    def addHydrogens(self):
        """
        fills valences with Hydrogens so that atom charge remains neutral
        """
        
        adjMatrix = self.adjMatrix
        atomNames = self.atomNames
        unpaireds = self.unpaireds
        lonepairs = self.lonepairs
        charges = self.charges
        for i, name in enumerate(self.atomNames):
            bondNum = sum(adjMatrix[i,:])
            if unpaireds:
                presentbds = bondNum+unpaireds[i]
            else:
                presentbds = bondNum
            bds = Adjacency.atomDict[name]
            bdif = bds-presentbds
            bdif = int(bdif)
            assert bdif >= 0, 'some atom has more bonds than it should'
            if bdif != 0:
                for k in xrange(bdif):
                    n = adjMatrix.shape[0]
                    adjMatrix = np.concatenate((adjMatrix.T,np.zeros((1,n)))).T
                    adjMatrix = np.concatenate((adjMatrix,np.zeros((1,n+1))))
                    adjMatrix[n,i] = 1
                    adjMatrix[i,n] = 1
                    atomNames.append('H')
                    if unpaireds:
                        unpaireds.append(0)
                    if lonepairs:
                        lonepairs.append(0)
                    if charges:
                        charges.append(0)
        
        self.adjMatrix = adjMatrix
        self.atomNames = atomNames
        self.unpaireds = unpaireds
        self.lonepairs = lonepairs
        self.charges = charges
            
        
    def getAdjList(self):
        """
        generates an RMG compatible adjacency list
        """
        
        adjlistStr = ''
        for i,name in enumerate(self.atomNames):
            newstr = str(i+1)+' '+name+' '
            if self.unpaireds:
                newstr += 'u'+str(int(self.unpaireds[i]))+ ' '
            else: #assume no radicals
                newstr += 'u0 '
            if self.lonepairs:
                newstr += 'p'+str(int(self.lonepairs[i]))+' '
            elif name != 'H':
                es = 8-2*sum(self.adjMatrix[i,:])
                if self.unpaireds:
                    es = es - self.unpaireds[i]
                lp = es // 2
                newstr += 'p'+str(int(lp))+' '
            elif name == 'H':
                newstr += 'p0 '
            if self.charges:
                newstr += 'c'+str(int(self.charges[i]))+' '
            for j in xrange(self.adjMatrix.shape[0]):
                if self.adjMatrix[i,j] != 0:
                    newstr += '{'+str(j+1)+','+str(self.adjMatrix[i,j])+'} '
            adjlistStr += newstr[:-1]+'\n'
            
        self.AdjList = adjlistStr
            
        return adjlistStr
        
    def toRMGmol(self):
        """
        generates an RMG molecule
        """
        self.getAdjList()
        self.RMGspc = Molecule().fromAdjacencyList(self.AdjList)
        return self.RMGspc

