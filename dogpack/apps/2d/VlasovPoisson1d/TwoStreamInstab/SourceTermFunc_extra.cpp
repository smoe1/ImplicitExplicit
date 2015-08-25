#include "dogdefs.h"

// This is a user-supplied routine that sets the
// source term at all the points "xpts"
void SourceTermFunc_extra(const dTensor2& xpts, const dTensor2& qvals, 
        const dTensor2& auxvals, dTensor2& source, void *data)
{
    int i,m;
    int numpts=xpts.getsize(1);
    int meqn=source.getsize(2);
    double x,v;
    double t = *(double*)data;
    double tmp

    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i,1);
        v = xpts.get(i,2);

        for (m=1; m<=meqn; m++)
        {
            tmp = 
                2.0*v*cos(2.0*x − 2.0*pi*t) − 2.0*pi*sin(2.0*x − 2.0*pi*t)  
                + 2.0*pi*v*sin(2.0*x − 2.0*pi*t)(2.0 − cos(2.0*x − 2.0*pi*t));
            tmp *= exp(-pow(2.0*v, 2) );
            source.set(i,m, tmp );
        }
    }

} 
