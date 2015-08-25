#include <cmath>
#include "dogdefs.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{

    const int numpts=xpts.getsize(1);
    const int meqn = 3;

    for(int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);

        const double r2 = pow(x-0.40,2)+pow(y-0.50,2);
        const double r  = sqrt(r2);

        for( int m=1; m <=meqn; m++ )
        {
            if (r<0.3)
            {  
                qvals.set(i,m, pow( cos(5.0/3.0*pi*r) ,6) );  
            }
            else
            {  qvals.set(i,m, 0.0 ); }
        }

        //qvals.set(i,1, sin(2*pi*x) );

    }
}
