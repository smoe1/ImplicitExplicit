///////////////////////////////////////////////////////////////////////////////
// Function that integrates the source term
///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"

///////////////////////////////////////////////////////////////////////////////
// This routine integrates psi.  The points provided in xpts are the quadrature
// points located at time t^{n+1}, and the time values provided in tpts are
// times in between time tn and tnew.
//
// This is the "backwards" integration of the source term.
///////////////////////////////////////////////////////////////////////////////
void IntegratePsiBackward( double tnew, const dTensor2 &xpts, const dTensor1 &tpts, dTensor3 &psi)
{

    const int numpts  = xpts.getsize(1);
    const int meqn    = psi.getsize(2);

    double a,x,y,v,t;

    for( int i=1; i <= numpts; i++ )
    for( int me=1; me <= meqn; me++ )
    {

        a = xpts.get(i,1);  // x-coordinate of quadrature point //
        v = xpts.get(i,2);  // v-coordinate of quadrature point //

        for( int n=1; n <= tpts.getsize(); n++ )
        {
            t   = tpts.get(n);
            x   = a+v*(t-tnew); // value along the characteristic at time t

            // source term from f_t + v f_x + E f_v = \psi //
            double ps = (1.0/2.0)*sin(-2.0*x+2.0*pi*t)*exp(-(1.0/4.0)*pow( (4.0*v-1.0),2))*(4.0*pi-4.0*v-8.0*sqrt(pi)*v+2.0*sqrt(pi)+4.0*sqrt(pi)*cos(-2.0*x+2.0*pi*t)*v-sqrt(pi)*cos(-2.0*x+2.0*pi*t));

            psi.set(i, me, n, ps);

        }
    }

}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// This routine supplies the "forwards" integration of the source term.
//
// This routine integrates psi.  The points provided in xpts are the quadrature
// points located at time t^{n}, and the time values provided in tpts are
// times in between time tn and tnew.
void IntegratePsiForward( double tn, const dTensor2 &xpts, const dTensor1 &tpts, dTensor3 &psi)
{

    const int numpts  = xpts.getsize(1);
    const int meqn    =  psi.getsize(2);

    double a,x,y,v,t;

    for( int i=1; i <= numpts; i++ )
    for( int me=1; me <= meqn; me++ )
    {

        a = xpts.get(i,1);  // x-coordinate of quadrature point //
        v = xpts.get(i,2);  // v-coordinate of quadrature point //

        for( int n=1; n <= tpts.getsize(); n++ )
        {

            t   = tpts.get(n);
            x   = a+v*(t-tn);

            // source term from f_t + v f_x + E f_v = \psi //
            double ps = (1.0/2.0)*sin(-2.0*x+2.0*pi*t)*exp(-(1.0/4.0)*pow( (4.0*v-1.0),2))*(4.0*pi-4.0*v-8.0*sqrt(pi)*v+2.0*sqrt(pi)+4.0*sqrt(pi)*cos(-2.0*x+2.0*pi*t)*v-sqrt(pi)*cos(-2.0*x+2.0*pi*t));

            psi.set(i, me, n, ps);

        }
    }

}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// This routine integrates psi.  The points provided in xpts are the quadrature
// points located at time t^{n+1}
//
// This is the "backwards" integration of the source term, and is performed
// exactly.
//////////////////////////////////////////////////////////////////////////////
void IntegratePsiBackwardExact( double tn, double dt, const dTensor2 &xpts, dTensor2 &psi)
{

    const int numpts  = xpts.getsize(1);
    const int meqn    = psi.getsize(2);

    for( int i=1; i <= numpts; i++ )
    for( int me=1; me <= meqn; me++ )
    {

        double x = xpts.get(i,1);  // x-coordinate of quadrature point //
        double v = xpts.get(i,2);  // v-coordinate of quadrature point //

//      double A1 = cos( 2.0*( x + v*dt - pi*tn) );
//      double A2 = cos( 2.0*(-x + pi*( tn+dt) ) );
//      double A3 = cos( 4.0*( x + v*dt - pi*tn ) );
//      double A4 = cos( 4.0*(-x + pi*( tn+dt) ) );

//      double M  = exp(-0.25* pow(4.0*v - 1.0, 2) ) / (16.0*(v+pi));
//      double c1 = -32.0*sqrt(pi)*v + 8.0* sqrt(pi) + 16.0*(pi-v);
//      double c2 = sqrt(pi) * ( 4.0*v-1.0 );

        double A1 = cos( 2.0*(-x + v*dt + pi*tn) );
        double A2 = cos( 2.0*(-x + pi*( tn+dt) ) );
        double A3 = cos( 4.0*(-x + v*dt + pi*tn ) );
        double A4 = cos( 4.0*(-x + pi*( tn+dt) ) );

        double M  = -exp(-0.25* pow(4.0*v - 1.0, 2) ) / (16.0*(-v+pi));
        double c1 = 32.0*sqrt(pi)*v - 8.0* sqrt(pi) - 16.0*(pi-v);
        double c2 = sqrt(pi) * ( 4.0*v-1.0 );


        double ps = M*(c1*(A1-A2) - c2*(A3-A4) );

        psi.set(i, me, ps);

    }

}
//////////////////////////////////////////////////////////////////////////////

