#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// source term at all the points "xpts"
void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals, 
        const dTensor2& auxvals, dTensor2& source)
{
    int i,m;
    int numpts=xpts.getsize(1);
    int meqn=source.getsize(2);
    double x,v, xmpit, tmp;

double t = 0.;

// TODO THIS NEEDS TO HAVE THE CORRECT TIME PLACED HERE!
printf("Error: the source term needs to have time placed correctly in here\n");
exit(1);

    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i,1);
        v = xpts.get(i,2);

        for (m=1; m<=meqn; m++)
        {
            xmpit = 2.0*(x-pi*t);
            tmp = ( 2.0*(v-pi)*sin(xmpit) 
                + 2.0*sqrt(pi)*v*sin(xmpit)*( 2.0-cos(xmpit)) ) * exp(-4.0*v*v);
            source.set(i,m, tmp );
        }
    }
}
