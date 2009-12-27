#!/usr/bin/python
#===============================================================================
#    Copyright (C) 2009
#    Babak Ahmadi [babak dot ahmadi at iais dot fraunhofer dot de]
#    Fabian Hadiji [fabian dot hadiji at iais dot fraunhofer dot de]
#    Kristian Kersting (coordination) [kristian dot kersting at iais dot fraunhofer dot de]
#
#    STREAM Project at
#        Fraunhofer IAIS, Sankt Augustin, Germany, and 
#        KDML, Unversity of Bonn, Germany 
#
#    This file is part of libSTREAM.
#
#    libSTREAM is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    libSTREAM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License 
#    along with this program; if not, see <http://www.gnu.org/licenses/>.
#===============================================================================
'''
Created on Dec 4, 2009

@author: Fabian Hadiji
'''
import sys;
import getopt
import math;
from libSTREAMWrapper import *;
import stream.Utils;

# add 'prior' factors first
def createInitialSubGraph (fg, remainingFactors):
    subGraphFactors = vector_factor();
    remainingVars = set();
    for var in fg.vars():
        remainingVars.add(var.label());
    
    for factorIdx, factor in enumerate(fg.factors()):
        if factor.vars().size() == 1:
            subGraphFactors.push_back(factor);
            remainingVars.remove(factor.vars().elements()[0].label());
            remainingFactors.remove(factorIdx);
    
    for factorIdx, factor in enumerate(fg.factors()):
        if factor.vars().size() > 1:
            add = False
            for var in factor.vars().elements():
                if var.label() in remainingVars:
                    remainingVars.remove(var.label());
                    add = True;
            if add:
                subGraphFactors.push_back(factor);
                remainingFactors.remove(factorIdx);
    
    assert fg.nrFactors() == (len(remainingFactors) + len(subGraphFactors))
    assert len(remainingVars) == 0;
    
    yield subGraphFactors;
    
    subGraph = FactorGraph(subGraphFactors);
    assert subGraph.nrVars() == fg.nrVars();
    
    yield subGraph;

# important: this calculation assumes:
# i) each factor is a conjunction of binary variables
# i) all factors have same weight (=1)
def calcGain(infAlg, factorIdx):
    # determine zero state:
    nonZeroState = -1
    zeroState = -1;
    factor = infAlg.factor(factorIdx);
    for i in range(factor.states()):
        if factor[i] == 1:
            zeroState = i;
        else:
            nonZeroState = i;
        if nonZeroState > -1 and zeroState > -1:
            break;
    
    
    states = calcState(factor.vars(), zeroState);
    p_f = 1;
    for pair in states:        
        varIdx = infAlg.findVar(pair.key())
        if isinstance(infAlg, CBP):
            belief = infAlg.beliefV(infAlg.reprV(varIdx));
        else:
            belief = infAlg.beliefV(varIdx)
        if pair.data() == 0:
            p_f *= belief[0];
        else:
            p_f *= belief[1];
    p_f = 1 - p_f;
    # weight calculation works for CNS:
    weight = math.log(factor[nonZeroState]);
    gain = math.log((1 - p_f) + p_f * math.exp(weight)) - p_f * weight;
    return gain;   

def runMarginalCpi(fgFilename, verbose):

    opts = PropertySet();
    opts.read('property_set_cbp');
    
    fg = FactorGraph();
    fg.ReadFromFile(fgFilename);
    
    bp = BP(fg, opts);
    bp.init();
    bp.run();
    
    # stores factor idx in the original graph
    remainingFactors = set();
    # stores factor idx in the original graph
    subgraphFactors = set();

    for idx in range(fg.nrFactors()):
        remainingFactors.add(idx);
       
    print "Original factor graph contains",fg.nrFactors(),"factors and",fg.nrVars(),"vars";
    
    factors, subGraph = createInitialSubGraph(fg, remainingFactors);

    if verbose > 0:
        print "initial sub-graph contains",subGraph.nrFactors(),"factors and",subGraph.nrVars(),"vars";
        
    infAlg = CBP(subGraph, opts);
    infAlg.init();
    infAlg.run();
    
    if verbose > 0:
        print "lifted sub-graph contains", infAlg.nrFactors(), "factors and",infAlg.nrVars(),"vars";
        print "maxnorm:",stream.Utils.maxnorm(bp, infAlg);

    eps = 0.01

    added = True;
    while added == True:
        added = False;
        subgraphFactors.clear()
        for idx in remainingFactors:
            gain = calcGain(bp, idx);
            if verbose > 2:
                print "gain:", gain
            if gain > eps:
                added = True;
                factors.push_back(fg.factor(idx))
                subgraphFactors.add(idx)
                
        remainingFactors -= subgraphFactors;
                          
        subgraph = FactorGraph(factors);
        infAlg = CBP(subgraph, opts);
        infAlg.init();
        infAlg.run();
        if verbose > 0:
            print "number of ground factors:",len(factors);
            print "sub-graph contains", infAlg.nrFactors(), "factors and",infAlg.nrVars(),"vars";
            print "maxnorm:",stream.Utils.maxnorm(bp, infAlg);
    
    if verbose > 1:
        print "marginals:"
        for i, var in enumerate(fg.vars()):
            print "  %-*d %-*s" % (5, var.label(), 35, infAlg.beliefV(i))

def usage():
    print "usage: runMarginalCpi.py [-v] [-h] [--help] fgFilename"
    print "-v"
    print "\tverbose level; default is 1"
    print "fgFilename"
    print "\tpath to fg file"

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hv:", ["help"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    
    verbose = 1;    
    for o, a in opts:
        if o == "-v":
            verbose = int(a)
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"
    
    if len(args) == 1:
        runMarginalCpi(args[0], verbose)
    else:
        usage()
        sys.exit(2)
        