#include "dogdefs.h"

// This is a user-supplied routine that sets the the boundary conditions
//
//      DOUBLE PERIODIC BOUNDARY CONDITIONS
//
void SetBndValues(dTensorBC4& q, dTensorBC4& aux)
{

    const int mx   = q.getsize(1);
    const int my   = q.getsize(2);
    const int meqn = q.getsize(3);
    const int kmax = q.getsize(4);
    const int mbc  = q.getmbc();
    const int maux = aux.getsize(3);

    int i,j,k,ell,m;
    double tmp;

    // -----------------------
    // BOUNDARY CONDITION LOOP
    // -----------------------
   
    for (ell=1; ell<=kmax; ell++)
    {
    
        // ***********************************************
        // LEFT BOUNDARY
        // ***********************************************
        for (i=0; i>=(1-mbc); i--)
        {
            for (j=1; j<=my; j++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(i+mx,j,m,ell);
                    q.set(i,j,m,ell, tmp );
                }
                
                // aux values
                for (m=1; m<=maux; m++)
                {
                    tmp = aux.get(i+mx,j,m,ell);
                    aux.set(i,j,m,ell, tmp );
                }
            }
        }
        // ***********************************************
        
        
        
        // ***********************************************
        // RIGHT BOUNDARY
        // ***********************************************
        for (i=(mx+1); i<=(mx+mbc); i++)
        {
            for (j=1; j<=my; j++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(i-mx,j,m,ell);
                    q.set(i,j,m,ell, tmp );
                }
                
                // aux values
                for (m=1; m<=maux; m++)
                {
                    tmp = aux.get(i-mx,j,m,ell);
                    aux.set(i,j,m,ell, tmp );
                }
            }
        }
        // ***********************************************
        
        
        
        // ***********************************************
        // BOTTOM BOUNDARY
        // ***********************************************
        for (j=0; j>=(1-mbc); j--)
        {
            for (i=1; i<=mx; i++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(i,j+my,m,ell);
                    q.set(i,j,m,ell, tmp );
                }               
                
                // aux values
                for (m=1; m<=maux; m++)
                {
                    tmp = aux.get(i,j+my,m,ell);
                    aux.set(i,j,m,ell, tmp );
                }
            }
        }
        // ***********************************************
        
        
        
        // ***********************************************
        // TOP BOUNDARY
        // ***********************************************
        for (j=(my+1); j<=(my+mbc); j++)
        {
            for (i=1; i<=mx; i++)
            {           
                // q values
                for (m=1; m<=meqn; m++)
                {
                    tmp = q.get(i,j-my,m,ell);
                    q.set(i,j,m,ell, tmp );
                }
                                
                // aux values
                for (m=1; m<=maux; m++)
                {
                    tmp = aux.get(i,j-my,m,ell);
                    aux.set(i,j,m,ell, tmp );
                }
            }
        }
        // ***********************************************
        
        
        // ***********************************************
        // BOTTOM LEFT CORNER
        // ***********************************************
        for (i=1; i<=mbc; i++)
          for (j=1; j<=mbc; j++)
            {
              for (m=1; m<=meqn; m++)
                {     
                  q.set(1-i,1-j,m,ell, q.get(mx+1-i,my+1-j,m,ell) );
                }
              for (m=1; m<=maux; m++)
                {     
                  aux.set(1-i,1-j,m,ell, aux.get(mx+1-i,my+1-j,m,ell) );
                }
            }
        // ***********************************************
        
        
        // ***********************************************
        // BOTTOM RIGHT CORNER
        // ***********************************************
        for (i=1; i<=mbc; i++)
          for (j=1; j<=mbc; j++)
            {
              for (m=1; m<=meqn; m++)
                {     
                  q.set(mx+i,1-j,m,ell, q.get(i,my+1-j,m,ell) );
                }
              for (m=1; m<=maux; m++)
                {     
                  aux.set(mx+i,1-j,m,ell, aux.get(i,my+1-j,m,ell) );
                }
            }
        // ***********************************************
        
        
        // ***********************************************
        // TOP RIGHT CORNER
        // ***********************************************
        for (i=1; i<=mbc; i++)
          for (j=1; j<=mbc; j++)
            {
              for (m=1; m<=meqn; m++)
                {     
                  q.set(mx+i,my+j,m,ell, q.get(i,j,m,ell) );
                }
              for (m=1; m<=maux; m++)
                {     
                  aux.set(mx+i,my+j,m,ell, aux.get(i,j,m,ell) );
                }
            }
        // ***********************************************
        
        
        // ***********************************************
        // TOP LEFT CORNER
        // ***********************************************
        for (i=1; i<=mbc; i++)
          for (j=1; j<=mbc; j++)
            {
              for (m=1; m<=meqn; m++)
                {     
                  q.set(1-i,my+j,m,ell, q.get(mx+1-i,j,m,ell) );
                }
              for (m=1; m<=maux; m++)
                {     
                  aux.set(1-i,my+j,m,ell, aux.get(mx+1-i,j,m,ell) );
                }
            }
        // ***********************************************
        
    }

}
