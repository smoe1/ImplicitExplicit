#include "tensors.h"

// Periodic boundary conditions
void SetBndValues(const dTensor2& node, 
        dTensorBC3& aux, 
        dTensorBC3& q)
{

    int i,m,ell;
    double tmp;
    int melems = q.getsize(1);
    int meqn   = q.getsize(2);
    int kmax   = q.getsize(3);
    int maux   = aux.getsize(2);
    int mbc    = q.getmbc();
    void L2Project(int,int,int,dTensor2,dTensorBC3,dTensorBC3,dTensorBC3&,
            void (*Func)(dTensor1, dTensor2, dTensor2, dTensor2&));


    // ***********************************************
    // LEFT BOUNDARY
    // ***********************************************
    for (i=0; i>=(1-mbc); i--)
    {        
        // q values
        for (m=1; m<=meqn; m++)
            for (ell=1; ell<=kmax; ell++)	  
            {
                tmp = q.get(i+melems,m,ell);
                q.set(i,m,ell, tmp );
            }

        // aux values
        for (m=1; m<=maux; m++)
            for (ell=1; ell<=kmax; ell++)	  
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
            for (ell=1; ell<=kmax; ell++)	  
            {
                tmp = q.get(i-melems,m,ell);
                q.set(i,m,ell, tmp );
            }

        // aux values
        for (m=1; m<=maux; m++)
            for (ell=1; ell<=kmax; ell++)
            {
                tmp = aux.get(i-melems,m,ell);
                aux.set(i,m,ell, tmp );
            }
    }
    // ***********************************************

}
