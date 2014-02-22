'''
Created on Feb 20, 2014

@author: Larry
'''
from importMol import importMol
husky=importMol("benz.mol")

print "Molecular Formula:"+husky.molecularFormula
print "charge:"+str(husky.molecularCharge)
print "Molar Mass:"+str(husky.molecularMass)+" g mol^-1"
