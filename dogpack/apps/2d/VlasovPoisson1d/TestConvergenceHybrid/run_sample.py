#!/usr/bin/env python
from __future__ import with_statement
from contextlib import closing
from subprocess import call, Popen, PIPE
import sys

dogpack_data_template = '''
; Parameters common to dogpack applications
[dogParams]
defaults_file = "$DOGPACK/config/dogParams_defaults.ini"
ndims       = 2          ; 1 or 2
mesh_type   = Cartesian  ; (either Cartesian or Unstructured) 
nout        = 1          ; number of output times to print results
tfinal      = 1.00       ; final time
dtv(1)      = %(dt)e      ; initial dt
dtv(2)      = %(dt)e      ; max allowable dt 
cflv(1)     = 5.10       ; max allowable Courant number
cflv(2)     = 5.00       ; desired Courant number
nv          = 500000     ; max number of time steps per call to DogSolve
time_stepping_method = User-Defined ; (e.g., Runge-Kutta, SDC, Lax-Wendroff, User-Defined)
;time_stepping_method = Runge-Kutta ; (e.g., Runge-Kutta, SDC, Lax-Wendroff, User-Defined)
limiter_method = moment ; (e.g., moment, viscosity)
space_order = %(s_order)i   ; =method(1)= order of accuracy in space
time_order  = %(t_order)i   ; =method(2)= order of accuracy in time
use_limiter = 0   ; =method(3)= use limiter (1-yes, 0-no)
verbosity   = 1   ; =method(4)= verbosity of output
mcapa       = 0   ; =method(5)= mcapa (capacity function index in aux arrays)
maux        = 2   ; =method(6)= maux (number of aux arrays, maux >= mcapa)
source_term = 1   ; =method(7)= source term (1-yes, 0-no)
meqn        = 1   ; number of equations
mrestart    = 0   ; restart from old data (1-yes, 0-no)
nstart      = 0   ; if mrestart==1: from file q(nstart)_restart.dat
datafmt     = 1   ; 1 for ascii, 5 for hdf5.
; -------------------------------------------------
;   Cartesian grid data (ignored if Unstructured)
; -------------------------------------------------
[grid]
mx    =  %(mx)i  ; number of grid elements in x-direction
my    =  %(my)i  ; number of grid elements in y-direction
mbc   = 2        ; number of ghost cells on each boundary

xlow  = -3.141592653589793
xhigh =  3.141592653589793
ylow  = -3.141592653589793
yhigh =  3.141592653589793

; -------------------------------------------------
;   Vlasov Parameters
; -------------------------------------------------
[vlasovparams]
;; Maxwellian Params: ;;
;;    fM(t=0,x,v) = rho_0 * exp(-0.5*v^2 / temp ) / sqrt(2*pi*temp) 
;;
;; Landau paramters:
;;
;;    f(t=0,x,v) = (1 + alpha*cos(kx) ) * fM(v)
;;
;alpha = 0.01
;k     = 0.5
'''

#defaults = {
#    'mx_ratio':3,
#    'iterations':10,
#    'ts_method':'Lax-Wendroff'
#}
#def main(mx_ratio=defaults['mx_ratio'],
#         iterations=default['iterations'],
#         ts_method=['ts_method']):
def main( iterations, mx_start, my_start, dt_start, n_start):
    '''runs time stepping method given by ts_method for 'iterations' number of
    times.  The resolution gets increased by a factor of 'mx_ratio' for each
    iteration.  The starting resolution is given by 'mx_start'.   Output is
    written to folders 'output0000n' for n=0...iterations-1.
'''
    data_file = 'parameters.ini'

    print 'Time order = (choose 2 or 4 )'
    s = sys.stdin.readline()
    time_order = int(s)

    print 'Space order = (choose 1,2,3,4 or 5 )'
    s = sys.stdin.readline()
    space_order = int(s)

    ratio = 2

    for i in range(iterations):
        mx_now = mx_start * ratio**i
        my_now = my_start * ratio**i
        dt_now = dt_start / ratio**i

        # we want to do:
        #   data = open('dogpack.data','w')
        #   print >> data, dogpack_data_template % { 'mx': mx_now, 'ts_method': ts_method} 
        #   data.close()
        # and we avoid the .close() (even in case of exception) with 'with':
        with closing(open(data_file,'w')) as data:
            # print >> data, dogpack_data_template % locals() 
            ## if we had used same names for everything
            print >> data, dogpack_data_template \
                % {'s_order':space_order, 't_order' : time_order, 'mx' : mx_now, 'my' : my_now, 'dt' : dt_now }
        # if you want to capture script output, do
        #   Popen(thing to run, shell=True, stdout=PIPE).communicate()[0]
        call('./dog.exe -o /var/tmp/seal/TestConvergenceHybrid/output%03i' % (i+n_start), shell=True)

if __name__ == '__main__':
    import optparse
    parser = optparse.OptionParser(
        usage='''%%prog (-h | [-m MX_RATIO] [-i ITERATIONS] [-x MX_START] [-w MY_START] 
                        [-t DT_START] [-n N_START])
    
%s''' % main.__doc__)

    parser.add_option('-i', '--iterations', type='int', default=4, 
                       help='''number of times we wish to run computation.
                       Default = 6''')
    parser.add_option('-x', '--mx-start', type='int', default=10, 
                       help='''MX_START = number of grid points to start
                       computation with''')
    parser.add_option('-w', '--my-start', type='int', default=10, 
                       help='''MY_START = number of grid points to start
                       computation with''')
    parser.add_option('-t', '--dt-start', type='float', default= 1.0,
                       help='''DT_START = starting dt value''')
    parser.add_option('-n', '--n-start', type='int', default=0, 
                       help='''N_START = folder number to start output from.
                       Folders are named output000i with 
                       i = N_START...N_START- 1.  The default value is 0''')
    opts, args = parser.parse_args()
#    check that the user provided enough arguments...
#    if not opts.infile or not opts.outfile:
#        parser.error('Both options -i and -o are required. Try -h for help.')
    main( opts.iterations, opts.mx_start,
        opts.my_start, opts.dt_start, opts.n_start)
