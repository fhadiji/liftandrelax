#!/usr/bin/python
#===============================================================================
#
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
Created on Nov 23, 2009

@author: Fabian Hadiji
'''

import os;
import sys;
from libSTREAMWrapper import *;

def maxState (belief):
    maxState = 0;
    maxVal = 0;
    for i, val in enumerate(belief):
        if val > maxVal:
            maxVal = val;
            maxState = i;
        
    return maxState;

def runCbpOnFgFile(fgFilename):
    print "Running on CBP on fg file", fgFilename
    opts = PropertySet();
    opts.read('property_set_cbp');

    print opts;

    cfg = CFactorGraph();
    cfg.readFromFile(fgFilename);
    
    print "cfg.nrVars:", cfg.nrVars();
    print "cfg.nrFactors:", cfg.nrFactors();
    print "cfg.nrEdges:", cfg.nrEdges();
    
    cbp = CBP(cfg,opts);
    cbp.init();
    
    print "cbp.nrVars:", cbp.nrVars();
    print "cbp.nrFactors:", cbp.nrFactors();
    print "cbp.nrEdges:", cbp.nrEdges();    
    
    cbp.run();        
    
    print "marginals:"
    for i, var in enumerate(cfg.vars()):
        print "  %-*d %-*s %d" % (5, var.label(), 35, cbp.beliefV(cbp.reprV(i)), maxState(cbp.beliefV(cbp.reprV(i))))


if __name__ == '__main__':
    if len(sys.argv) == 2:
        runCbpOnFgFile(sys.argv[1])
    else:
        print "usage: runCbpOnFgFile.py <fgFilename>"