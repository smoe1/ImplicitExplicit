#!/usr/bin/env python
import numpy as np
import sys #used for command line arguments

def check_positivity( outputdir ):

    outputdir = outputdir + '/'

    try:
        fid = open( outputdir + 'qhelp.dat', 'r' )
        mesh_type   = fid.readline()
        meqn        = int(   fid.readline() )
        maux        = int(   fid.readline() )
        nout        = int(   fid.readline() )
        space_order = int(   fid.readline() )
        mx          = int(   fid.readline() )
        my          = int(   fid.readline() )
        xlow        = float( fid.readline() )
        xhigh       = float( fid.readline() )
        xlow        = float( fid.readline() )
        xhigh       = float( fid.readline() )
        fid.close()
    except IOError as e:
        print 'ERROR - could not open folder: ' + outputdir.strip('/')
        return
    
    'derived quantities from qout'
    kmax = (int) ( space_order*(space_order+1)/2 )

    'check the positivity of each output file provided'
    for n in range(0, nout+1):
        output_num = n

        filename  = 'q%04d' % output_num + '.dat'
        data = np.loadtxt( outputdir + '/' + filename )

        time = data[0]
        data = data[1:]

        print 'checking %d points provided in file: ' % (mx*my) \
            + filename
        q1   = data[:mx*my]
        for n in range(0, len(q1) ):
            if( q1[n] < -1e-10 ):
                print 'you went negative!  Q1 = %2.8e' % q1[n]

if __name__ == '__main__':

    if( len(sys.argv) < 2 ):
        outputdir = 'output'
    else:
        outputdir = sys.argv[1]
    check_positivity( outputdir )
