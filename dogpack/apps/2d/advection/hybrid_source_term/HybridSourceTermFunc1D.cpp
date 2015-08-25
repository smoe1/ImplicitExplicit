#include "dogdefs.h"
#include <cmath>

// This is a user-supplied routine that sets the source term
// This is a dummy function, used for hybrid stepping that should be overwritten
// by examples wishing to do so
//
// TODO - this v really shouldn't be hard coded here ...
//
void HybridSourceTermFunc1D(double t, const dTensor1& speeds, const dTensor1& xpts, dTensor2& fvals)
{
    int numpts = xpts.getsize();
    int meqn   = fvals.getsize(2);

//printf("  t = %f\n", t );

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);
        for (int me=1; me<=meqn; me++)
        {
            double v =  speeds.get(me);
            //fvals.set(i,me, 0.5*cos(2*pi*x) );
            //fvals.set(i,me, 0.5*cos(2*pi*v) );
            //fvals.set(i,me, 0.5 + 0.5*cos(2*pi*t) );

            double tmp = -pi*cos(x)*cos(pi*t)-pi*sin(x)*sin(pi*t)+v*cos(x)*cos(pi*t)+v*sin(x)*sin(pi*t);
            fvals.set(i,me,tmp);
        }
    }        

}
