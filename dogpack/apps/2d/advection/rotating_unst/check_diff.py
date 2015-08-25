#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import os
from math import log
import numpy as np

def main( dirs ):
    '''Compute the difference between scalars stored in two different files.
'''

    folder_parallel = (os.getcwd() + '/' + dirs[0] )
    folder_serial   = (os.getcwd() + '/' + dirs[1] )

    # TODO - check for the final frame, or perhaps be nice and ask the user 
    # which frame number she wants to use.
    q_parallel  = np.loadtxt(folder_parallel + "/q0010.dat")[1:]
    q_serial    = np.loadtxt(folder_serial   + "/q0010.dat")[1:]
    
    # Compute the l2-norm of all the coefficients:
    diff = q_parallel-q_serial
    diff = np.sqrt( np.dot(diff,diff) )

    print(diff)

if(__name__ == '__main__'):

    import argparse

    # User help message
    parser = argparse.ArgumentParser(description='Check the difference between two runs')

    # names for the two directories
    parser.add_argument('--out1', default='output_parallel', help='First  directory')
    parser.add_argument('--out2', default='output_serial',   help='Second directory')

    # parse the arguments:
    args = parser.parse_args()

    # check the l2-norm of these two directories
    main( [args.out1, args.out2] )
