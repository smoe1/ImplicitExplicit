#include "dogdefs.h"
#include "DogParams.h"
#include "DogParamsCart2.h"

#include <cmath>
#include "VlasovParams.h"   // for maxwellian parameters

void L2Project(int istart, int iend, int jstart, int jend,
        const dTensorBC4& q,
        const dTensorBC4& aux, dTensorBC4& Fout,
        void (*Func)(const dTensor2& xpts,
            const dTensor2& qvals,
            const dTensor2& auxvals,
            dTensor2& source));

void RHSBoundaryFunc(const dTensor2& xpts,
            const dTensor2& qvals,
            const dTensor2& auxvals,
            dTensor2& fvals)
{

    const int numpts=xpts.getsize(1);
    const double rho = vlasovParams.rho0;
    const double T   = vlasovParams.temp;

    for (int i=1; i<=numpts; i++)
    {
        double x = xpts.get(i,1);
        double v = xpts.get(i,2);
        double tmp = rho / sqrt( 2.0*pi*T ) * exp( -0.5 * v*v / T );
        fvals.set(i, 1,  tmp);       
    }


}

// This is a user-supplied routine that sets the the boundary conditions
//
//      The scheme assumes zero outflow boundary conditions.  Here, one uses
//      the function RHSBoundaryFunc to define what the right hand side boundary
//      condition should be.  It is assumed that f_bndy = f_bndy( y ) only.
//
//      ZEROTH ORDER EXTRAPOLATION IN Y DIRECTION
//
//      NOTHING IS SET FOR THE FOUR CORNERS OF THE DOMAIN
//
void SetBndValues(dTensorBC4& q, dTensorBC4& aux)
{

    const int mx   = q.getsize(1);  assert_eq( mx, dogParamsCart2.get_mx());
    const int my   = q.getsize(2);  assert_eq( my, dogParamsCart2.get_my());
    const int meqn = q.getsize(3);  assert_eq( meqn, dogParams.get_meqn());
    const int kmax = q.getsize(4);  assert_eq( kmax, dogParams.get_kmax());
    const int mbc  = q.getmbc();    assert_eq( mbc, dogParamsCart2.get_mbc());
    const int maux = aux.getsize(3);

    dTensorBC4 rhs(1, my, meqn, kmax, mbc );
    L2Project(1, 1, 1, my, q, aux, rhs, &RHSBoundaryFunc );


    // -----------------------
    // BOUNDARY CONDITION LOOP
    // -----------------------
    for (int ell=1; ell<=kmax; ell++)
    {      

        // ***********************************************
        // BOTTOM BOUNDARY (Periodic)
        // ***********************************************
        for (int j=0; j>=(1-mbc); j--)
            for (int i=1; i<=mx; i++)
            {           
                // q values
                for (int m=1; m<=meqn; m++)
                {
                    double tmp = q.get(i,-j+1,m,ell);
                    q.set(i,j,m,ell, tmp );
                }               

                // aux values
                for (int m=1; m<=maux; m++)
                {
                    double tmp = aux.get(i,-j+1,m,ell);
                    aux.set(i,j,m,ell, tmp );
                }
            }
        // ***********************************************      

        // ***********************************************
        // TOP BOUNDARY (Periodic)
        // ***********************************************
        for (int j=(my+1); j<=(my+mbc); j++)
            for (int i=1; i<=mx; i++)
            {           
                // q values
                for (int m=1; m<=meqn; m++)
                {
                    double tmp = q.get(i, 2*my-j+1,m,ell);
                    q.set(i,j,m,ell, tmp );
                }

                // aux values
                for (int m=1; m<=maux; m++)
                {
                    double tmp = aux.get(i, 2*my-j+1,m,ell);
                    aux.set(i,j,m,ell, tmp );
                }            
            }
        // ***********************************************

        // ***********************************************
        // LEFT BOUNDARY (Dirchlet)
        // ***********************************************
        for (int i=0; i>=(1-mbc); i--)
        for (int j=1; j<=my; j++)
        {           
            // q values
            for (int m=1; m<=meqn; m++)
            {
                q.set(i,j,m,ell, 0. );
            }

            // aux values
            for (int m=1; m<=maux; m++)
            {
                double tmp = aux.get(i+mx,j,m,ell);
                aux.set(i,j,m,ell, tmp );
            }
        }
        // ***********************************************


        // ***********************************************
        // RIGHT BOUNDARY (Dirchlet)
        // ***********************************************
        for (int i=(mx+1); i<=(mx+mbc); i++)
        for (int j=1; j<=my; j++)
        {           
            // q values
            for (int m=1; m<=meqn; m++)
            {
                q.set(i,j,m,ell, rhs.get(1,j,m,ell) );
            }

            // aux values
            for (int m=1; m<=maux; m++)
            {
                double tmp = aux.get(i-mx,j,m,ell);
                aux.set(i,j,m,ell, tmp );
            }
        }
        // ***********************************************



        // ***********************************************
        // BOTTOM LEFT CORNER
        // ***********************************************
//      for (int i=1; i<=mbc; i++)
//          for (int j=1; j<=mbc; j++)
//          {
//              for (int m=1; m<=meqn; m++)
//              {     
//                  q.set(1-i,1-j,m,ell, q.get(mx+1-i,my+1-j,m,ell) );
//              }
//              for (int m=1; m<=maux; m++)
//              {     
//                  aux.set(1-i,1-j,m,ell, aux.get(mx+1-i,my+1-j,m,ell) );
//              }
//          }
        // ***********************************************


        // ***********************************************
        // BOTTOM RIGHT CORNER
        // ***********************************************
//      for (int i=1; i<=mbc; i++)
//          for (int j=1; j<=mbc; j++)
//          {
//              for (int m=1; m<=meqn; m++)
//              {     
//                  q.set(mx+i,1-j,m,ell, q.get(i,my+1-j,m,ell) );
//              }
//              for (int m=1; m<=maux; m++)
//              {     
//                  aux.set(mx+i,1-j,m,ell, aux.get(i,my+1-j,m,ell) );
//              }
//          }
        // ***********************************************


        // ***********************************************
        // TOP RIGHT CORNER
        // ***********************************************
//      for (int i=1; i<=mbc; i++)
//          for (int j=1; j<=mbc; j++)
//          {
//              for (int m=1; m<=meqn; m++)
//              {     
//                  q.set(mx+i,my+j,m,ell, q.get(i,j,m,ell) );
//              }
//              for (int m=1; m<=maux; m++)
//              {     
//                  aux.set(mx+i,my+j,m,ell, aux.get(i,j,m,ell) );
//              }
//          }
        // ***********************************************


        // ***********************************************
        // TOP LEFT CORNER
        // ***********************************************
//      for (int i=1; i<=mbc; i++)
//          for (int j=1; j<=mbc; j++)
//          {
//              for (int m=1; m<=meqn; m++)
//              {     
//                  q.set(1-i,my+j,m,ell, q.get(mx+1-i,j,m,ell) );
//              }
//              for (int m=1; m<=maux; m++)
//              {     
//                  aux.set(1-i,my+j,m,ell, aux.get(mx+1-i,j,m,ell) );
//              }
//          }
        // ***********************************************


    }

}
