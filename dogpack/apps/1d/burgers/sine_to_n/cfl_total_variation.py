def main( ):
    '''Open folders output00[nums] and extract the total variation and conserved
    quatities from these files.

    See also: ConSoln.cpp.
'''

    import numpy as np

    # parser for the parameters.ini files (needed to read in the cfl number)
    import ConfigParser, os
    config = ConfigParser.ConfigParser()

    # Plotting flag
    plt_soln = 1

    # Count the number of output directories:
    my_dictionary = {}
    i = 1
    while( 1 ):

        directory_num = my_dictionary['dir_num'] =  i

        folder = (os.getcwd() + '/output%(dir_num)04i/') % my_dictionary
        if( not os.path.exists(folder) ):
            break
        i = i+1

    tvd_vec = np.zeros( (i,1) )
    cfl_vec = np.zeros( (i,1) )

    print("found %i folders" % len(tvd_vec) )

    print("-----------------------------------------------------------------------")
    print("CFL; Deviation in TV; Largest undershoot; Largest overshoot")
    print("-----------------------------------------------------------------------")

    for i in range( len(tvd_vec) ):

        directory_num = my_dictionary['dir_num'] = i
        folder = (os.getcwd() + '/output%(dir_num)04i/') % my_dictionary

        # Grab the current cfl number:
        config.readfp(open(folder + '/parameters.ini'))
        cfl_now = float( config.get('dogParams','cflv(2)') )

        data = np.loadtxt( folder + '/total-variation.dat' )
#       data = np.loadtxt( folder + '/conservation.dat' )
        t    = data[:,0]
        tv   = data[:,1]
        minq = min( data[:,2] )
        maxq = max( data[:,3] )
        dtv  = max( tv[1:] - tv[0:len(tv)-1] )
        print("%2.3f;    %2.5e;    %2.5e;    %2.5e" % ( cfl_now, dtv, min(minq,0.), max(maxq-1.0,0.) ) )

        tvd_vec[i] = dtv
        cfl_vec[i] = cfl_now

    if( plt_soln ):

        import numpy as np
        import matplotlib.pyplot as plt

        fig = plt.figure()

        ax1  = fig.add_subplot(111)
        #ax1.plot( cfl_vec, tvd_vec, 'bo' )
        ax1.semilogy( cfl_vec, tvd_vec, 'bo' )
        ax1.set_xlim( (min(cfl_vec), max( cfl_vec ) ) )
        ax1.set_ylim( (-0.1*max( abs(tvd_vec) ), 1.1*max( tvd_vec )) ) 

        ax1.set_xlabel('CFL number')
        ax1.set_ylabel('Increase in Total Variation')

        #ax2  = fig.add_subplot(221)
        #ax2.plot( cfl_good, tvd_good, 'bo' )
        #ax2.plot( cfl_bad, tvd_bad, 'kx' )
        #ax2.set_xlim( (min(cfl_good), max( cfl_good ) ) )
        #ax2.set_ylim( (-0.1, max( max( tvd_good ), max( tvd_bad ) )+0.1 ) )

        fig.show()

if __name__ == '__main__':
    main.__doc__

    main( )
