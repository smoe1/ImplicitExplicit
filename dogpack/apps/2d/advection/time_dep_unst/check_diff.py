#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import os
from math import log
import numpy as np

def main( ):
    '''Write some help documentation here
'''

    folder_parallel = (os.getcwd() + '/output_parallel/')
    folder_serial   = (os.getcwd() + '/output_serial/')

    q_parallel  = np.loadtxt(folder_parallel + "/q0010.dat")[1:]
    q_serial    = np.loadtxt(folder_serial   + "/q0010.dat")[1:]
      
    diff = q_parallel-q_serial

    diff = np.sqrt( np.dot(diff,diff) )

    print diff

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser(
        usage='''%%prog (-h |
    
%s''' % main.__doc__)

    opts, args = parser.parse_args()

    main( )
