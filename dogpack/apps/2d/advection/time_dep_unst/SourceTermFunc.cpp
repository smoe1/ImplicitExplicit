#include "dogdefs.h"
#include <cmath>

// This is a user-supplied routine that sets the
// source term at all the points "xpts"
void SourceTermFunc(const dTensor2& xpts, const dTensor2& qvals, 
        const dTensor2& auxvals, dTensor2& source)
{

    const double t = qvals.get(1,2);

    const int numpts=xpts.getsize(1);
    const int meqn=source.getsize(2);
    for (int i=1; i<=numpts; i++)
    {

        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);

        // Center of circle:
        const double x0 = 0.0;
        const double y0 = 0.0;

        // each component of the velocity:
        const double vx = auxvals.get(i,1);
        const double vy = auxvals.get(i,2);

        double r  = sqrt( pow(x-x0,2) + pow(y-y0,2) );

        double dr_dx = (x-x0) / (r+1e-13);
        double dr_dy = (y-y0) / (r+1e-13);

        // ------------------------ //
        // Exact solution is given by:
        // qex = (0.75 + 0.25*cos(2*pi*t) ) * ( 0.5 + 0.5*cos( pi*r ) )
        // ------------------------ //
        double qt = -0.5*pi*sin(2.0*pi*t)   * ( 0.5 + 0.5*cos( pi*r ) );
        double qx = (0.75+0.25*cos(2*pi*t)) * ( -0.5*pi*sin( pi*r )* dr_dx );
        double qy = (0.75+0.25*cos(2*pi*t)) * ( -0.5*pi*sin( pi*r )* dr_dy );
        
        source.set(i,1, qt + vx*qx + vy*qy );

        // eqn = 2 is used to store time
        source.set(i,2, 1.0e0 );
    }
}
