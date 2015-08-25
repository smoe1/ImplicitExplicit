#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import os
from math import log

run_matlab_template = """
addpath('%(dogpack_dir)s/matlab/');
run_convergence2( %(num_pts)i, '%(curr_folder)s' );
"""

def main( ):
    '''Write some help documentation here
'''
    input_file = 'tmpfile.m'
    dogpack_home_dir = os.environ['DOGPACK']
    number_pts = 1
    iterations = 4

    my_dictionary = {'dogpack_dir' : dogpack_home_dir, 'num_pts' : number_pts}
    print "# leading comments can be given a '#' character"
    old_err = i = 0
    while( 1 ):

        directory_num = my_dictionary['dir_num'] =  i
#       folder = 'outputSL%(s_order)i_%(t_order)i_00%(dir_num)i' % my_dictionary
        folder = '/var/tmp/seal/TestConvergence/output00%(dir_num)i' % my_dictionary
        if( not os.path.exists(folder) ):
            break

        my_dictionary['curr_folder'] = folder
        # we want to do:
        #   data = open('dogpack.data','w')
        #   print >> data, dogpack_data_template % { 'mx': mx_now, 'ts_method': ts_method} 
        #   data.close()
        # and we avoid the .close() (even in case of exception) with 'with':
        directory_num = i
        with closing(open(input_file,'w')) as data:
            # print >> data, dogpack_data_template % locals() 
            ## if we had used same names for everything
            print >> data, run_matlab_template % my_dictionary
       
        s = Popen('matlab -nodesktop -nosplash < ' + input_file, shell=True, stdout=PIPE).communicate()[0]
        new_err = float( s[s.find('>>')-1:len(s)].replace(">>","") )

        r1 = '%(new).3e  ' % {'old': old_err, 'new' : new_err}
        
        if( old_err > 0 and new_err > 0 ):
            result = r1 + folder + '   log2(ratio) = %(rat).3f' % \
                {'rat' : log( (old_err/new_err), 2) } 
        else:
            result = r1 + folder + '   log2(ratio) = %(rat).3f' % \
                {'old' : old_err, 'new' : new_err, 'rat' : (old_err/new_err) } 
        
        print result

        old_err = new_err

        i = i + 1

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser(
        usage='''%%prog (-h |
    
%s''' % main.__doc__)

    opts, args = parser.parse_args()

    main( )
