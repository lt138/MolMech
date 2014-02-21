'''
Created on Feb 20, 2014

@author: Larry
'''
import numpy as np, itertools as its, math as m

class importMol():
    
    def __init__(self, molFile):
        
        self.readMolFile(molFile)
        self.bondLength(self.cartMatrix,self.bondMatrix,self.numBond,self.numAtom)
        self.bondAngle(self.cartMatrix,self.bondMatrix,self.numBond,self.numAtom,self.bondLengthM)
    
    def readMolFile(self,molFile):
        """reads an input mdl mol file line by line, determining the number of bonds and atoms in the molecule
        
        ,a matrix of Cartesian coordinates for each atom, and a matrix with the bond connectivity for each atom"""
        
        mF=open(molFile, 'r')
        lines=mF.readlines()
        #3rd line of mol file contains number of atoms and number of bonds
        atomCounts = lines[3].split()
        self.numAtom = int(atomCounts[0])
        self.numBond = int(atomCounts[1])
        
        self.atomType=[]
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
                    self.bondLengthM[x,y]=np.linalg.norm(np.subtract(cartMatrix[x,:],cartMatrix[y,:]))
                    
    def bondAngle(self,cartMatrix,bondMatrix,numBond,numAtom,bondLengthM):
        """ Determines the bond angle between all the groups of 3 bonded atoms in the molecule"""
       
        # enumerate all the possible ways to arrange atoms in the molecule into groups of 3 without repetition. 
        comb= list(its.combinations(range(numAtom),3))
        self.bondAngleM=np.zeros((len(comb),4))
        
        for x in range(0,len(comb)):
            #only compute the bond angle for bonded groups of 3 atoms.
            if (bondMatrix[comb[x][1],comb[x][0]]!= 0) and (bondMatrix[comb[x][1],comb[x][2]]!= 0):
                #determine the distance vectors along the internuclear axis and takes the dot product of the two 
                vec1=np.subtract(cartMatrix[comb[x][0],:],cartMatrix[comb[x][1],:])
                vec2=np.subtract(cartMatrix[comb[x][2],:],cartMatrix[comb[x][1],:])
                dotp=np.dot(vec1,vec2)
                # computes the bond angle by taking the arccosine of the dot product of the two vectors
                # divided by the product of the vector norms
                self.bondAngleM[x,0]=m.acos(dotp/(bondLengthM[comb[x][1],comb[x][0]]*bondLengthM[comb[x][1],comb[x][2]]))
                self.bondAngleM[x,1:]=[comb[x][0],comb[x][1],comb[x][2]]
                
            else:
                self.bondAngleM[x,1:]=[comb[x][0],comb[x][1],comb[x][2]] 