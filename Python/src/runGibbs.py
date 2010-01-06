#!/usr/bin/python

import os;
import sys;
from libSTREAMWrapper import *;

def runGibbs(fgFilename):
    opts = PropertySet();
    opts.read('property_set_cbp');

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
    
    opts.Set("iters","5000")
    gibbs = Gibbs(fg,opts);
    #gibbs.init();
    gibbs.run();   

    opts.Set("iters","100000")
    exact = Gibbs(fg,opts);
    #gibbs.init();
    exact.run();   


    #jtreeOpts = PropertySet();
    #jtreeOpts.Set("updates","HUGIN")
    #jtreeOpts.Set("verbose","1")

    #jtree = JTree(fg,jtreeOpts);
    #jtree.init();
    #jtree.run();
    
    print "marginals:"

    avgError = 0
    avgLiftedError = 0
    for lifted in range(cbp.nrVars()):
        print "Lifted:", lifted
        avg = 0.0
        for ground in range(len(cbp.clusterV(lifted))):
            orig = cbp.clusterV(lifted)[ground]
            belief = gibbs.beliefV(orig)
            avg = avg + belief[0]
        avg = avg / len(cbp.clusterV(lifted)) 
        for ground in range(len(cbp.clusterV(lifted))):
            orig = cbp.clusterV(lifted)[ground]
            belief = gibbs.beliefV(orig)
            exactBelief = exact.beliefV(orig)
            error = abs(belief[0]-exactBelief[0])
            liftedError = abs(avg - exactBelief[0])
            print "  %-*d %-*s %-*s %-*s %-*s %-*s" % (5, ground,5, 1, 35, belief, 35, exactBelief, 5, error, 5, liftedError)
            avgError = avgError + error
            avgLiftedError = avgLiftedError + liftedError
        print "Avg:", avg
            #print "Ground:", cbp.clusterV(lifted)[ground]
    print "Avg. Error:", avgError / fg.nrVars()
    print "Avg. Lifted Error:", avgLiftedError / fg.nrVars()

    #print gibbs.beliefV(0)
    #for i, var in enumerate(fg.vars()):
    #    print "  %-*d %-*s %-*s" % (5, var.label(),5, 1, 35, gibbs.beliefV(var.label()))


if __name__ == '__main__':
    if len(sys.argv) == 2:
        runGibbs(sys.argv[1])
    else:
        print "usage: runGibbs.py <fgFilename>"
