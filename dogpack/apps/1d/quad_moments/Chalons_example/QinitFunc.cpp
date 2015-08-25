#include "dogdefs.h"
#include <cmath>
#include "QuadMomentParams.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
// Each application is REQUIRED to define one of these.
//
// Input:
//
//    xpts( 1:numpts )           - The x-coordinates for a list of points
//
// Output:
//
//    qvals( 1:numpts, 1:meqn )  - The vector of conserved variables, q at
//                                 each point.
//
// See also: AuxFunc.
void QinitFunc(const dTensor1& xpts, dTensor2& qvals)
{

    const int numpts = xpts.getsize();

    // Initial conditions
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);
        double rho,u,p,q,r;


        // equilibrium distribution with constant density       

        /*
           rho = 1.0;
           u   = 0.0;
           p   = rho;
           q   = 0.0;
           r   = 3.0*rho;
         */

        // equilibrium distribution with rho = sqrt(2*pi)/2*(2+cos(2*pi*x))
        /*
           rho = sqrt(2.0*pi)/2.0*(2.0+cos(2*pi*x));
           u   = 0.0;
           p   = rho;
           q   = 0.0;
           r   = 3.0*rho;
         */  


        if(x<0.5e0)
        {
            rho =  1.0;
            u   =  1.0;
            p   =  1.0/3.0;
            q   =  0.0;
            r   =  0.3;//2.7*p*p/rho;
        }
        else 
        { 
            rho =  1.0;
            u   = -1.0;
            p   =  1.0/3.0;
            q   =  0.0;
            r   =  0.3;//2.7*p*p/rho;
        }

        /*
           if(x<0.5e0)
           {
           rho =  1.0;
           u   =  0.0;
           p   =  1.0;
           q   =  0.0;
           r   =  3.0*p*p/rho;
           }
           else 
           { 
           rho =  0.125;
           u   =  0.0;
           p   =  0.1;
           q   =  0.0;
           r   =  3.0*p*p/rho;
           }
         */


        // Check input
        const double p2 = p*p;
        const double p3 = p*p2;
        const double q2 = q*q;
        const double rmin = (p2/rho + q2/p);
        if (rho<=0.0 || p<=0.0 || r<rmin )
        {
            printf("\n");   
            printf(" ERROR in Initial Conditions: \n");
            printf("  rho = %22.16e \n",rho);
            printf("    p = %22.16e \n",p);
            printf("    q = %22.16e \n",q);
            printf("    r = %22.16e \n",r);
            printf("\n");
            printf(" rmin = %22.16e \n",rmin);
            printf("\n");
            exit(1);
        }

        const double u2 = u*u;
        const double u3 = u2*u;
        const double u4 = u3*u;

        const double M2 = p + rho*u2;
        const double M3 = q + 3.0*p*u + rho*u3;
        const double M4 = r + 4.0*q*u + 6.0*p*u2 + rho*u4;

        qvals.set(i,1, rho   );  // density
        qvals.set(i,2, rho*u );  // momentum
        qvals.set(i,3, M2    );  // energy
        qvals.set(i,4, M3    );  // 3rd-order moment
        qvals.set(i,5, M4    );  // 4th-order moment 
    }

}
