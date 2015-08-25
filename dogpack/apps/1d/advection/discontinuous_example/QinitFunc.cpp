#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"

// parameters for "4 bells" problem
const double a = 0.5;
const double z = -0.7;
const double delta = 0.005;
const double alpha = 10.0;
const double beta  = log(2.0)/(36.*delta*delta);

inline double G( double xi, double xc )
{   return exp(-beta*(xi-xc)*(xi-xc)); }

inline double F( double xi, double xc )
{   return sqrt( Max( 1.0-alpha*alpha*(xi-xc)*(xi-xc), 0. ) ); }

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
    const int numpts=xpts.getsize();

    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);

        double tmp = 0.;
        if(-0.8 <= x && x <= -0.6 )
        {
            tmp = ( G(x, z-delta) + G(x,z+delta) + 4.0*G(x,z) ) / 6.0;
        }        
        else if( -0.4 <= x && x <= -0.2 )
        { tmp = 1.0; }
        else if( 0. <= x && x <= 0.2 )
        {
            tmp = 1.0 - fabs( 10.*(x-0.1) );
        }
        else if( 0.4 <= x && x <= 0.6 )
        {
            tmp = ( F(x, a-delta) + F(x,a+delta) + 4.0*F(x,a) ) / 6.0;
        }
        qvals.set(i,1, tmp );
    }

}
