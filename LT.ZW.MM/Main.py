'''
Created on Feb 20, 2014

@author: Larry
'''
from importMol import importMol
husky=importMol("meth.mol")
print husky.numAtom
print husky.numBond
print husky.atomType
print husky.bondMatrix
print husky.bondLengthM
print husky.bondAngleM
