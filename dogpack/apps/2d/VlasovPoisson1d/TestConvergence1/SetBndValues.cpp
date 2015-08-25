#include "dogdefs.h"

// This is a user-supplied routine that sets the the boundary conditions
//
//      PERIODIC BOUNDARY CONDITIONS IN X DIRECTION
//
//      ZEROTH ORDER EXTRAPOLATION IN Y DIRECTION
//
//      NOTHING IS SET FOR THE FOUR CORNERS OF THE DOMAIN
//
void SetBndValues(dTensorBC4& q, dTensorBC4& aux)
{
    int mx   = q.getsize(1);
    int my   = q.getsize(2);
    int meqn = q.getsize(3);
    int kmax = q.getsize(4);
    int mbc  = q.getmbc();
    int maux = aux.getsize(3);

    // -----------------------
    // BOUNDARY CONDITION LOOP
    // -----------------------
    for (int ell=1; ell<=kmax; ell++)
    {      
        // ***********************************************
        // BOTTOM BOUNDARY
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
        // TOP BOUNDARY
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

        // ***********************************************
        // LEFT BOUNDARY (Periodic)
        // ***********************************************
        for (int i=0; i>=(1-mbc); i--)
            for (int j=1; j<=my; j++)
            {           
                // q values
                for (int m=1; m<=meqn; m++)
                {
                    double tmp = q.get(i+mx,j,m,ell);
                    q.set(i,j,m,ell, tmp );
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
        // RIGHT BOUNDARY (Periodic)
        // ***********************************************
        for (int i=(mx+1); i<=(mx+mbc); i++)
            for (int j=1; j<=my; j++)
            {           
                // q values
                for (int m=1; m<=meqn; m++)
                {
                    double tmp = q.get(i-mx,j,m,ell);
                    q.set(i,j,m,ell, tmp );
                }

                // aux values
                for (int m=1; m<=maux; m++)
                {
                    double tmp = aux.get(i-mx,j,m,ell);
                    aux.set(i,j,m,ell, tmp );
                }
            }
        // ***********************************************


    }

}
