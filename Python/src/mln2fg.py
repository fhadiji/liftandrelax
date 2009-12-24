#!/usr/bin/python
#===============================================================================
#    Copyright (C) 2009  
#    Babak Ahmadi [babak dot ahmadi at iais dot fraunhofer dot de]
#    Fabian Hadiji [fabian dot hadiji at iais dot fraunhofer dot de]
#    Supervisor: Kristian Kersting [kristian dot kersting at iais dot fraunhofer dot de]
# 
#    Fraunhofer IAIS, STREAM Project, Sankt Augustin, Germany
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
Created on Oct 23, 2009

@author: Fabian Hadiji
'''
import sys
import getopt
from os import remove
from libSTREAMWrapper import *
from MLN import *
from stream.Utils import generateFilename;
from stream.cbp.MlnToFg import *

def usage():
    print "usage: mln2fg.py [--cnf] [-o fgFilename] mlnFilename dbFilename"

def mln2fg(mlnFilename, dbFilename, outFilename, cnf):
    
    mlnToFg = MlnToFg();
       
    fg = mlnToFg.mlnFileToFactorGraph(mlnFilename, dbFilename, cnf);
    
    if outFilename == None:
        suffix = ".fg"
        if cnf:
            suffix = "_cnf.fg"
        outFilename = generateFilename(mlnFilename, suffix)
            
    open(outFilename,'w').write(fg.__str__())
    print 'Wrote factorgraph to', outFilename;
    

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:", ["help", "output=", "cnf"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
        
    cnf = False
    outFilename = None;
    for o, a in opts:
        
        if o == "--cnf":
            cnf = True;
        elif o in ("-o", "--output"):
            outFilename = a
        elif o in ("-h", "--help"):
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"
    
    if len(args) == 1:
        mln2fg(args[0], None, outFilename, cnf)
    elif len(args) == 2:
        mln2fg(args[0], args[1], outFilename, cnf)
    else:
        usage()
        sys.exit(2)

