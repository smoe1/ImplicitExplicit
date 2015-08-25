#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE

dogpack_data_template = '''
; Parameters common to dogpack applications
[dogParams]
defaults_file = "$DOGPACK/config/dogParams_defaults.ini"
ndims       = 1             ; 1 or 2
mesh_type   = Cartesian     ; (either Cartesian or Unstructured) 
nout        = 1             ; number of output times to print results
tfinal      = 1.0           ; final time
dtv(1)      = 1.0           ; initial dt
dtv(2)      = 1e10          ; max allowable dt 
cflv(1)     = 1.0           ; max allowable Courant number
cflv(2)     = %(cfl)f       ; desired Courant number
nv          = 500000     ; max number of time steps per call to DogSolve
time_stepping_method = %(ts_method_str)s ; (e.g., Runge-Kutta, SDC, Lax-Wendroff, User-Defined)
limiter_method = moment ; (e.g., moment, viscosity)
space_order = %(s_order)i   ; =method(1)= order of accuracy in space
time_order  = %(t_order)i   ; =method(2)= order of accuracy in time
use_limiter = 0   ; =method(3)= use limiter (1-yes, 0-no)
verbosity   = 1   ; =method(4)= verbosity of output
mcapa       = 0   ; =method(5)= mcapa (capacity function index in aux arrays)
maux        = 1   ; =method(6)= maux (number of aux arrays, maux >= mcapa)
source_term = 0   ; =method(7)= source term (1-yes, 0-no)
meqn        = 1   ; number of equations
mrestart    = 0   ; restart from old data (1-yes, 0-no)
nstart      = 0   ; if mrestart==1: from file q(nstart)_restart.dat
datafmt     = 1   ; 1 for ascii, 5 for hdf5.
withPyClawPlotting = 0; (1-yes, 0-no)

[grid]
mx    =  %(mx)i  ; number of grid elements in x-direction
mbc   =      2   ; number of ghost cells on each boundary
xlow  =  0.0e0   ; left end point
xhigh =  1.0e0   ; right end point
'''

def main(cfl, ts_method, space_order, time_order, iterations, mx_start, n_start):
    '''Simple script for performing a batch refinement study.

  This script runs time stepping method given by ts_method for 'iterations' number of
  times.  The resolution gets increased by a factor of 'mx_ratio' for each
  iteration.  The starting resolution is given by 'mx_start'.   Output is
  written to folders 'output0000n' for n=0...iterations-1.
'''
    data_file = 'parameters.ini'
    ratio = 2

    integrators   = ['Runge-Kutta', 'SDC', 'Lax-Wendroff', 'User-Defined']
    ts_method_str = integrators[ts_method]
    print(ts_method_str)

    for i in range(iterations):
        mx_now = mx_start * ratio**i

        # we want to do:
        #   data = open('dogpack.data','w')
        #   print >> data, dogpack_data_template % { 'mx': mx_now, 'ts_method': ts_method} 
        #   data.close()
        # and we avoid the .close() (even in case of exception) with 'with':
        with closing(open(data_file,'w')) as data:

            # print >> data, dogpack_data_template % locals() 
            ## if we had used same names for everything

            my_dictionary = {'s_order' : space_order, 't_order' : time_order, 
                    'cfl' : cfl,
                    'mx' : mx_now,
                    "i_now": (i+n_start),
                    'ts_method_str' : ts_method_str }
            print >> data, dogpack_data_template % my_dictionary

        # if you want to capture script output, do
        #   Popen(thing to run, shell=True, stdout=PIPE).communicate()[0]
        #%cmd = './dog.exe -o outputSL%(s_order)i_%(t_order)i_%(i_now)03i' % my_dictionary
        cmd = './dog.exe -o output00%(i_now)i' % my_dictionary
        print cmd
        call(cmd, shell=True)
        print ''' 
//////////////////////// finished running output directory output%03i //////////
''' % (i+n_start)

def parse_input( help_message ):
  
    import argparse, sys
  
    parser = argparse.ArgumentParser (
      prog = 'python '+ sys.argv[0],
      description =    help_message,
      formatter_class = argparse.RawTextHelpFormatter,
    )

    parser.add_argument('CFL',
                      type = float,
                      help = 'maximum Courant number in domain')
  
    parser.add_argument('-s','--time_integrator',
                      type    = int,
                      choices = range(4),
                      default =  0,
                      dest    = 't_stepper',
                      metavar = 'X',
                      help    = 
    '''choose integrator X for time-stepping:
  0. Runge-Kutta
  1. SDC
  2. Lax-Wendroff
  3. 'User-Defined' time integrator.  (See Makefile for what gets linked to)
(default: 0)''')

    parser.add_argument('-f','--frames',
                      type    = int,
                      nargs   = 2,
                      default = [10, 7],
                      metavar = ('MX_START', 'N_FRAMES'),
                      help    = 
''' Refinment parameters:
Produce N_FRAMES of refinements, starting with 
MX_START grid points.
(default: [10, 7] )''')

    parser.add_argument('-t', '--order',
                      type = int,
                      nargs   = 2,
                      default = [4, 4],
                      metavar = ('S_ORDER', 'T_ORDER'),
                      help = 
'''Order of accuracy in space S_ORDER, and time T_ORDER.
(default: [4,4])''')
  
    return parser.parse_args()

if __name__ == '__main__':
    """ Batch refinement study. 
    
    This script runs multiple samples, which can be later analyzed for a
    convergence study."""

    # Parse input arguments
    args = parse_input( main.__doc__ )
    print(args)
    print('')

    main( args.CFL, args.t_stepper, 
        args.order[0], args.order[1], args.frames[1], args.frames[0], 0 )
