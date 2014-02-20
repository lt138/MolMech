'''
Created on Feb 20, 2014

@author: Larry
'''
import numpy as np
import itertools as its
import math as m
class importMol():
    numAtom=0
    numBond=0
    cartMatrix=np.matrix(0)
    bondMatrix=np.matrix(0)
    atomType=[]
    bondLength=np.matrix(0)
    bondAngle=np.matrix(0)
    
    def __init__(self, molFile):
        self.readMolFile(molFile)
        self.bondLength(self.cartMatrix,self.bondMatrix,self.numBond,self.numAtom)
        self.bondAngle(self.cartMatrix,self.bondMatrix,self.numBond,self.numAtom,self.bondLengthM)
    
    def readMolFile(self,molFile):
        """reads an input mdl mol file line by line, determining the number of bonds and atoms in the molecule
        
        ,a matrix of Cartesian coordinates for each atom, and a matrix with the bond connectivity for each atom"""

        mF=open(molFile, 'r')
        lines=mF.readlines()
        
        atomCounts = lines[3].split()
        self.numAtom = int(atomCounts[0])
        self.numBond = int(atomCounts[1])
        
        self.cartMatrix  = np.zeros((self.numAtom,3))
        self.bondMatrix  = np.zeros((self.numAtom,self.numAtom))
        
        for x in range(4,4+self.numAtom):
            atomrow=lines[x].split()
            for y in range(0,2):
                self.cartMatrix[x-4,y]=atomrow[y]
            self.atomType.append(atomrow[3])
            
        for x in range(4+self.numAtom,4+self.numAtom+self.numBond):
            bondrow=lines[x].split()
            i=int(bondrow[0])-1
            j=int(bondrow[1])-1
            self.bondMatrix[i,j]= int(bondrow[2])
            self.bondMatrix[j,i]= int(bondrow[2]) 
            
        mF.close()
        
    def bondLength(self,cartMatrix,bondMatrix,numBond,numAtom):
        """Determines the between each atom using the Cartesian distance formula
        Takes the Cartesian coordinate matrix, bonding matrix, number of bonds and atoms as required inputs."""
        self.bondLengthM =np.zeros((numAtom,numAtom))
        for x in range(0,numAtom):
            for y in range(0,numAtom):
                if bondMatrix[x,y]!=0:
                    self.bondLengthM[x,y]=((cartMatrix[x,0]-cartMatrix[y,0])**2+(cartMatrix[x,1]-cartMatrix[y,1])**2+(cartMatrix[x,2]-cartMatrix[y,2])**2)**0.5
    def bondAngle(self,cartMatrix,bondMatrix,numBond,numAtom,bondLengthM):
        comb= list(its.combinations(range(numAtom),3))
        self.bondAngleM=np.zeros((len(comb),4))
        for x in range(0,len(comb)):
            if (bondMatrix[comb[x][1],comb[x][0]]!= 0) & (bondMatrix[comb[x][1],comb[x][2]]!= 0):
                vec1=np.subtract(cartMatrix[comb[x][1],:],cartMatrix[comb[x][0],:])
                vec2=np.subtract(cartMatrix[comb[x][1],:],cartMatrix[comb[x][2],:])
                dotp=np.dot(vec1,vec2)
                self.bondAngleM[x,0]=m.acos(dotp/(bondLengthM[comb[x][1],comb[x][0]]*bondLengthM[comb[x][1],comb[x][2]]))
                self.bondAngleM[x,1:]=[comb[x][0],comb[x][1],comb[x][2]]
            else:
                self.bondAngleM[x,1:]=[comb[x][0],comb[x][1],comb[x][2]] 