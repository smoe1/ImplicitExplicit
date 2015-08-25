#include <iostream>
#include <math.h>
#include "dogdefs.h"
#include "DogParams.h"

// This is a user-supplied routine that sets the
// auxiliary arrays at all the points "xpts"
//
//

double fdisc(double x,
	     double y)
{
    double fdisc;

    if(x>0.0 && y<0.55*x)
    {fdisc=1.0;}
    else{fdisc=-1.0;}
    return fdisc;
}


void AuxFunc(const dTensor2& xpts, 
	     dTensor2& auxvals)
{

    const int numpts = xpts.getsize(1);

    const int mpoints1d = dogParams.get_space_order();

    // Using space_order is safer than expecting sqrt to return the correct
    // integer.
//  dTensor1 w1d(pow(numpts,0.5));
//  dTensor1 x1d(pow(numpts,0.5));
    dTensor1 w1d(mpoints1d);
    dTensor1 x1d(mpoints1d);

    double cl     = 1.0;
    double rhol   = 1.0;
    double cr     = 0.5;
    double rhor   = 1.0;
    double bulkl  = cl*cl * rhol;
    double bulkr  = cr*cr * rhor;

    // Extract the Gaussian quadrature weights
    void setGaussPoints1d(dTensor1& w1d, dTensor1& x1d);
    setGaussPoints1d( w1d, x1d );

    // rho averages
    double cav      = 0.0;
    double rhoav    = 0.0;
    for (int i=1; i<=mpoints1d; i++)
    {
        for(int j=1;j<=mpoints1d;j++)
        {
            int   j1  = i+(mpoints1d-1)*j;
            double f  = fdisc( xpts.get(j1,1), xpts.get(j1,2) );
            cav      += 0.25*w1d.get(i)*w1d.get(j)*( 0.5*(f+1.0)*cr   + 0.5*(1.0-f)*cl   );
            rhoav    += 0.25*w1d.get(i)*w1d.get(j)*( 0.5*(f+1.0)*rhor + 0.5*(1.0-f)*rhol );
        }
    }

    for (int i=1; i<=numpts; i++)
    {  
        auxvals.set(i,1, cav   );
        auxvals.set(i,2, rhoav );
    }

}
