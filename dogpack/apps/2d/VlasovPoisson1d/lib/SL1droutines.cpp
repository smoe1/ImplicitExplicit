///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "dogdefs.h"
#include "dog_math.h"
#include "DogParams.h"
#include "DogParamsCart2.h"

// Periodic boundary conditions
void SetPeriodicBndValues1D( dTensorBC3& q)
{
    int i,m,ell;
    double tmp;
    const int melems = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);
    const int mbc    = q.getmbc();
    
    for (ell=1; ell<=kmax; ell++)
    { 
        // ***********************************************
        // LEFT BOUNDARY
        // ***********************************************
        for (i=0; i>=(1-mbc); i--)
        {        
	        // q values
	        for (m=1; m<=meqn; m++)
	        {
	            tmp = q.get(i+melems,m,ell);
		        q.set(i,m,ell, tmp );
            }
                
        }
        // ***********************************************  

        // ***********************************************
        // RIGHT BOUNDARY
        // ***********************************************
        for (i=(melems+1); i<=(melems+mbc); i++)
        {        
	        // q values
	        for (m=1; m<=meqn; m++)
	        {
	            tmp = q.get(i-melems,m,ell);    
		        q.set(i,m,ell, tmp );
            }
                
        }
        // ***********************************************
    }
}

//////////////////////////////////////////////////////////////////////////////
// Compute q_x using finite differences of moments. 
// Output is written to qout.
//
//    qout(:,1,:) = q_x
//    qout(:,2,:) = q_xx
//
//////////////////////////////////////////////////////////////////////////////
void FiniteDiff(int mestart, int meend, const dTensorBC3& q, dTensorBC3& qout)
{
    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);
    const int mbc    = q.getmbc();

    const double dx     = dogParamsCart2.get_dx();
    const double two_dx = 2.0*dx;
    const double dx2    = dx*dx;

    for (int i=1; i<=mx; i++)
    for( int me=mestart; me <= meend; me=me+2 )
    {
        dTensor1 dA ( 5  );
        dTensor1 dA2( 5  );
        dA.setall   ( 0. );
        dA2.setall  ( 0. );


        int ip1 = i+1;
        int im1 = i-1;
        int ic  = i;

        if (i<2)
        {
            im1 = im1+mx;
        }
        if (i>(mx-1))
        {
            ip1 = ip1-mx;
        }

        for (int k=1; k<=kmax; k++)
        {
            dA.set(  k, (q.get(ip1,1,k)-q.get(im1,1,k))/two_dx );
            dA2.set( k, (q.get(ip1,1,k)-2.0*q.get(ic,1,k)+q.get(im1,1,k))/dx2 );
        }

        switch(kmax)
        {
            case 5:
                qout.set(i,me,  5,  dA.get(5) );
                qout.set(i,me+1,5, dA2.get(5) );

            case 4:
                qout.set(i,me,  4,  dA.get(4) );
                qout.set(i,me+1,4, dA2.get(4) );

            case 3:
                qout.set(i,me,  3,  dA.get(3) - 14.0*sq5* dA.get(5) );
                qout.set(i,me+1,3, dA2.get(3) - 7.0*sq5*dA2.get(5) );

            case 2:
                qout.set(i,me,  2,  dA.get(2) - 10.0*sq7/sq3* dA.get(4) );
                qout.set(i,me+1,2, dA2.get(2) - 5.0*sq7/sq3*dA2.get(4) );

            case 1:
                qout.set(i,me,  1,  dA.get(1) - 2.0*sq5* dA.get(3) + 78.0* dA.get(5) );
                qout.set(i,me+1,1, dA2.get(1) - sq5*dA2.get(3) + 11.0*dA2.get(5) );

                break;
        }

    }

}
