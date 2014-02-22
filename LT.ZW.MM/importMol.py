'''
Created on Feb 20, 2014

@author: Larry
'''
import numpy as np, itertools as its, math as m,atomicParam as ap

class importMol():
    
    def __init__(self, molFile):
        
        self.readMolFile(molFile)
        self.getBondLength(self.cartMatrix,self.bondMatrix,self.numBond,
                           self.numAtom)
        self.getBondAngle(self.cartMatrix,self.bondMatrix,self.numBond,
                          self.numAtom,self.bondLength)
        self.getMolecularFormula(self.atomType)
        self.getMolecularMass(self.atomType,self.isotopeDiff)
        self.getMolecularCharge(self.atomCharge)
    
    def readMolFile(self,molFile):
        """reads an input mdl mol file line by line, determining the number
         of bonds and atoms in the molecule, a matrix of Cartesian
        coordinates for each atom, and a matrix with the bond
        connectivity for each atom"""
        
        mF=open(molFile, 'r')
        lines=mF.readlines()
        #3rd line of mol file contains number of atoms and number of bonds
        atomCounts=lines[3].split()
        self.numAtom=int(atomCounts[0])
        self.numBond=int(atomCounts[1])
        #A list containing the elements of the molecules as ordered
        #in the mol file.
        self.atomType=[]
        #Lists any deviations in element atomic mass listed in the mol file.
        self.isotopeDiff=[]
        #Lists any charges on an atom in the same order as atomType
        self.atomCharge=[]
        self.cartMatrix=np.zeros((self.numAtom,3))
        self.bondMatrix=np.zeros((self.numAtom,self.numAtom))
        
        for x in range(4,4+self.numAtom):
            
            atomrow=lines[x].split()
            
            for y in range(0,2):
                self.cartMatrix[x-4,y]=atomrow[y]
                
            self.atomType.append(atomrow[3])
            self.isotopeDiff.append(int(atomrow[4]))
            self.atomCharge.append(int(atomrow[5]))
            
        for x in range(4+self.numAtom,4+self.numAtom+self.numBond):
            
            bondrow=lines[x].split()
            
            i=int(bondrow[0])-1
            j=int(bondrow[1])-1
            
            self.bondMatrix[i,j]=int(bondrow[2])
            self.bondMatrix[j,i]=int(bondrow[2]) 
            
        mF.close()
        
    def getBondLength(self,cartMatrix,bondMatrix,numBond,numAtom):
        """Determines the between each atom using the Cartesian distance 
        formula. Takes the Cartesian coordinate matrix, bonding matrix,
        number of bonds and atoms as required inputs."""
        
        self.bondLength =np.zeros((numAtom,numAtom))
        
        for x in range(0,numAtom):
            for y in range(0,numAtom):
                if bondMatrix[x,y]!=0:
                    self.bondLength[x,y]=(np.linalg.norm(np.subtract(
                                    cartMatrix[x,:],cartMatrix[y,:])))
                    
    def getBondAngle(self,cartMatrix,bondMatrix,numBond,numAtom,bondLength):
        """ Determines the bond angle between all the groups of 3 bonded
         atoms in the molecule"""
       
        # enumerate all the possible ways to arrange atoms in
        # the molecule into groups of 3 without repetition. 
        comb= list(its.combinations(range(numAtom),3))
        self.bondAngle=np.zeros((len(comb),4))
        
        for x in range(0,len(comb)):
            #only compute the bond angle for bonded groups of 3 atoms.
            if ((bondMatrix[comb[x][1],comb[x][0]]!= 0) and
                   (bondMatrix[comb[x][1],comb[x][2]]!= 0)):
                #determine the distance vectors along the internuclear axis
                #and takes the dot product of the two 
                vec1=(np.subtract(cartMatrix[comb[x][0],:],
                                 cartMatrix[comb[x][1],:]))
                vec2=(np.subtract(cartMatrix[comb[x][2],:],
                                  cartMatrix[comb[x][1],:]))
                dotp=np.dot(vec1,vec2)
                # computes the bond angle by taking the arccosine
                # of the dot product of the two vectors divided by
                # the product of the vector norms
                self.bondAngle[x,0]=(m.acos(dotp/(bondLength[comb[x][1],
                                    comb[x][0]]*bondLength[comb[x][1],
                                                          comb[x][2]])))
                self.bondAngle[x,1:]=[comb[x][0],comb[x][1],comb[x][2]]
                
            else:
                self.bondAngle[x,1:]=[comb[x][0],comb[x][1],comb[x][2]] 
                
    def getMolecularFormula(self,atomType):
        self.molecularFormula=""
        
        #sorts the atomType list alphabetically without duplicate elements
        sortedAtomType=sorted(list(set(atomType)))

        #append element and count (if not equal to 1) to
        #the molecular formula in alphabetical order
        for x in range(0,len(sortedAtomType)):
            self.molecularFormula+=sortedAtomType[x]
            if atomType.count(sortedAtomType[x])!=1:
                self.molecularFormula+=str(atomType.count(sortedAtomType[x]))
                
    def getMolecularCharge(self,atomCharge):
        self.molecularCharge=sum(atomCharge)
        
    def getMolecularMass(self,atomType,isotopeDiff):
        self.molecularMass=0
        for x in range(0,len(atomType)):
            self.molecularMass+=isotopeDiff[x]+ap.getAtomicMass(atomType[x])