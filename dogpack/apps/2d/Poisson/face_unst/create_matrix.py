#!/usr/bin/env python

# template for running the unstructured code
#
# This should be consistent with whatever is currenlty located in
# Unstructured_Mesh/input2D.data.
#
# Replace as necessary
#
input2d_template = '''
Unstructured          : Grid type (Structured vs. Unstructured)
%(h0)f                : Initial grid spacing (H0)
50000                 : Maximum allowed iterations
%(sorder)d       : Create a submesh with given integer refinement factor
----------------------------------------------------------------
Bounding box:
-1.0e0                : xmin
 1.0e0                : xmax
-1.0e0                : ymin
 1.0e0                : ymax
----------------------------------------------------------------
4                    : Number of fixed points
----------------------------------------------------------------
Fixed points (x y):
    -1.000000000000000e+00    -1.000000000000000e+00
    -1.000000000000000e+00     1.000000000000000e+00
     1.000000000000000e+00    -1.000000000000000e+00
     1.000000000000000e+00     1.000000000000000e+00
'''

# matlab script template
matlab_script_template = '''
addpath('../../lib/matlab');
CreateMatrixGlobal( %(sorder)d, 0 );
exit;
'''

def main(sorder, h0):

    from contextlib import closing
    from subprocess import call, Popen, PIPE

    my_dictionary = {'h0':h0, 'sorder':sorder}

    # Create new mesh:
    with closing(open('Unstructured_Mesh/input2D.data','w')) as data:
        print >> data, input2d_template % my_dictionary
    cmd = 'cd Unstructured_Mesh; make; ./mesh.exe'
    call(cmd, shell=True)

    ## Create script file to call matlab with:
    # Create new mesh:
    with closing(open('matlab/script.m','w')) as data:
        print >> data, matlab_script_template % my_dictionary
    cmd = 'cd matlab; matlab -nodesktop -nosplash -nojvm -r script'
    call(cmd, shell=True)

if __name__=='__main__':

    import argparse
    parser = argparse.ArgumentParser( description='''Create a new mesh, call the
    matalab code responsible for setting up matrices for Poisson solver.''')

    parser.add_argument('-s', '--space_order', type=int, default=3, \
        help="Spatial order of accuracy.  (default=3)" )
    parser.add_argument('-h0', '--h0', type=float, default=0.029, \
        help="Initial grid spacing for unstructured mesh.  (default=0.029)" )

    args = parser.parse_args()

    main( args.space_order, args.h0 )
