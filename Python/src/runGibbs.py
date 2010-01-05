#!/usr/bin/python

import os;
import sys;
from libSTREAMWrapper import *;

def runGibbs(fgFilename):
    opts = PropertySet();
    opts.read('property_set_cbp');
    opts.Set("iters","100")

    print opts;

    fg = FactorGraph();
    fg.ReadFromFile(fgFilename);

    print "cfg.vars:", fg.vars();
    print "cfg.nrVars:", fg.nrVars();
    print "cfg.nrFactors:", fg.nrFactors();
    print "cfg.nrEdges:", fg.nrEdges();
        

    #created compressed network
    cfg = CFactorGraph(fg)
    cbp = CBP(cfg,opts);
    cbp.init();
    
    print "cbp.nrVars:", cbp.nrVars();
    print "cbp.nrFactors:", cbp.nrFactors();
    print "cbp.nrEdges:", cbp.nrEdges();    
    

    gibbs = Gibbs(fg,opts);
    #gibbs.init();
    gibbs.run();        
    
    print "marginals:"

    for lifted in range(cbp.nrVars()):
        print "Lifted:", lifted
        for ground in range(len(cbp.clusterV(lifted))):
            print "  %-*d %-*s %-*s" % (5, ground,5, 1, 35, gibbs.beliefV(cbp.clusterV(lifted)[ground]))
            #print "Ground:", cbp.clusterV(lifted)[ground]

    #print gibbs.beliefV(0)
    #for i, var in enumerate(fg.vars()):
    #    print "  %-*d %-*s %-*s" % (5, var.label(),5, 1, 35, gibbs.beliefV(var.label()))


if __name__ == '__main__':
    if len(sys.argv) == 2:
        runGibbs(sys.argv[1])
    else:
        print "usage: runGibbs.py <fgFilename>"
