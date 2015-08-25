#include "dogdefs.h"
#include <cmath>

void HybridSourceTermFunc1D(double t, const dTensor1& speeds, const dTensor1& xpts, dTensor2& fvals)
{
    int numpts = xpts.getsize();
    int meqn   = fvals.getsize(2);

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i);
        for (int me=1; me<=meqn; me++)
        {

            double v =  speeds.get(me);

            double tmp = 0.5*sin(2.*(x-pi*t))*exp(-0.25*pow(4.*v-1.,2) ); 

            tmp *= ( (2.*sqrt(pi)+1.)*(4.*v-2.*sqrt(pi)) -
                sqrt(pi)*(4.*v-1.)*cos(2.*(x-pi*t)));

            fvals.set(i,me, tmp);

        }
    }        

}
