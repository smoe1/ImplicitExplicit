#----------------------------------------------------------
def read_q1d_data(qfile):

    import string
    import numpy as np

    # open file
    Rqfile = open(qfile,'r')

    # get time
    linestring = Rqfile.readline()
    linelist = string.split(linestring)
    time = float(linelist[0])

    # get number of grid points
    linestring = Rqfile.readline()
    linelist = string.split(linestring)
    numpts = int(linelist[0])

    # get number of equations
    linestring = Rqfile.readline()
    linelist = string.split(linestring)
    meqn = int(linelist[0])
    
    # get mesh information
    r1d = np.zeros(numpts,float);
    for k in range (0,numpts):
        linestring = Rqfile.readline()
        linelist = string.split(linestring)
        r1d[k] = float(linelist[0])

    # get solution information
    q1d = np.zeros((numpts,meqn),float);
    for m in range (0,meqn):
        for k in range (0,numpts):
            linestring = Rqfile.readline()
            linelist = string.split(linestring)
            q1d[k,m] = float(linelist[0])
    
    # close file
    Rqfile.close()

    # return time
    return time,r1d,q1d
#----------------------------------------------------------
