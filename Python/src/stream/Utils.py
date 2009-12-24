'''
Created on Dec 4, 2009

@author: stream
'''
from libSTREAMWrapper import *;

def maxnorm(infAlg1, infAlg2):
    maxVal = -1;
    for i in range(infAlg1.nrVars()):
        for j in range(infAlg1.beliefV(i).states()):
            if isinstance(infAlg1, CBP):
                belief1 = infAlg1.beliefV(infAlg1.reprV(i));
            else:
                belief1 = infAlg1.beliefV(i);
                
            if isinstance(infAlg2, CBP):
                belief2 = infAlg2.beliefV(infAlg2.reprV(i));
            else:
                belief2 = infAlg2.beliefV(i);
            val = abs(belief1[j] - belief2[j])
            if val > maxVal:
                maxVal = val;
    
    return maxVal

def generateFilename(filename, newSuffix):
    newFilename = filename + newSuffix;
    rpos = filename.rfind('.')
    lpos = filename.rfind('/'); 
    if rpos > 0:
        newFilename = filename[lpos + 1:rpos] + newSuffix;    
    return newFilename;

def getMaxState (belief):
    maxState = -1;
    maxVal = -1;
    for i, val in enumerate(belief):
        if val > maxVal:
            maxState = i;
            maxVal = val;
    return maxState

def checkSolution (solution, fg):
    for factor in fg.factors():
        mapping = map_var_sizet();
        for var in factor.vars():
            mapping[var] =  solution[var.label()]
        state = factor.vars().calcState(mapping);
        # assuming exponentiated factors
        if factor[state] == 1:
            return False;

    return True;