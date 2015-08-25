#include "dogdefs.h"

// Periodic boundary conditions
void SetBndValues(const dTensor2& node, 
		  dTensorBC3& aux, 
		  dTensorBC3& q)
{
    int i,m,ell;
    double tmp;
    const int melems = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);
    const int maux   = aux.getsize(2);
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
                
	        // aux values
	        for (m=1; m<=maux; m++)
            {
	            tmp = aux.get(i+melems,m,ell);
		        aux.set(i,m,ell, tmp );
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
                
	        // aux values
	        for (m=1; m<=maux; m++)
            {
	            tmp = aux.get(i-melems,m,ell);
		        aux.set(i,m,ell, tmp );
            }
		}
        // ***********************************************
    }
}

/*
// Solid wall boundary conditions
void SetBndValues(const dTensor2& node, 
        dTensorBC3& aux, 
        dTensorBC3& q)
{
    const int melems = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);
    const int maux   = aux.getsize(2);
    const int mbc    = q.getmbc();

    for (int ell=1; ell<=kmax; ell++)
    { 
        // ***********************************************
        // LEFT BOUNDARY
        // ***********************************************
        for (int i=0; i>=(1-mbc); i--)
        {        
            // q values
            for (int m=1; m<=meqn; m++)
            {
                double tmp = q.get(1-i,m,ell);                    
                q.set(i,m,ell, tmp );
            }
            // normal velocity to wall
            double tmp = -q.get(1-i,2,ell);
            q.set(i,2,ell, tmp );

            // aux values
            for (int m=1; m<=maux; m++)
            {
                double tmp = aux.get(1,m,ell);                    
                aux.set(i,m,ell, tmp );
            }
        }
        // ***********************************************  


        // ***********************************************
        // RIGHT BOUNDARY
        // ***********************************************
        for (int i=(melems+1); i<=(melems+mbc); i++)
        {        
            // q values
            for (int m=1; m<=meqn; m++)
            {
                double tmp = q.get(2*melems-i,m,ell);                    
                q.set(i,m,ell, tmp );
            }
            // normal velocity to wall
            double tmp = -q.get(2*melems-i,2,ell);
            q.set(i,2,ell, tmp );

            // aux values
            for (int m=1; m<=maux; m++)
            {
                double tmp = aux.get(melems,m,ell);                    
                aux.set(i,m,ell, tmp );
            }
        }
        // ***********************************************

    }

}
*/
