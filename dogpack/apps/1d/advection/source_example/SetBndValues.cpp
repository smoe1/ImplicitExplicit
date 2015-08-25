#include<cmath>
#include "dogdefs.h"

// Boundary conditions
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
    void L2Project(int mopt, int istart, int iend,
            const dTensor2& node,
            const dTensorBC3& qin, 
            const dTensorBC3& auxin,  
            dTensorBC3& Fout,
            void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));
    void QBndLeft(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void AuxBndLeft(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);


    // ***********************************************
    // LEFT BOUNDARY --  Dirichlet boundary condition
    // ***********************************************
    L2Project(0,1-mbc,0,node,q,aux,q,&QBndLeft);      // q-value
    L2Project(0,1-mbc,0,node,q,aux,aux,&AuxBndLeft);  // aux-value
    // ***********************************************  


    // ***********************************************
    // RIGHT BOUNDARY -- OUTFLOW BOUNDARY
    // ***********************************************
    for (i=(melems+1); i<=(melems+mbc); i++)
    { 
        for (ell=1; ell<=kmax; ell++)
        {
            // q values
            for (m=1; m<=meqn; m++)
            {
                tmp = q.get(i-1,m,ell);

                q.set(i,m,ell, tmp );
            }

            // aux values
            for (m=1; m<=maux; m++)
            {
                tmp = aux.get(i-1,m,ell);

                aux.set(i,m,ell, tmp );
            }        	    
        }
    }
    // ***********************************************    
}


// ***********************************************
// Left Boundary - Q-values
void QBndLeft(const dTensor1& xpts,
        const dTensor2& NOT_USED_1, 
        const dTensor2& NOT_USED_2, 
        dTensor2& qout)
{
    int i;
    int numpts=xpts.getsize();
    double x;

    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i);

        qout.set(i,1, 1.0 + sin(2.0*pi*x) );
        //qout.set(i,1, exp(x) );
    }
}
// ***********************************************


// ***********************************************
// Left Boundary - Aux-values
void AuxBndLeft(const dTensor1& xpts, 
        const dTensor2& NOT_USED_1, 
        const dTensor2& NOT_USED_2, 
        dTensor2& auxout)
{
    int i;
    int numpts=xpts.getsize();
    double x;

    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i);

        auxout.set(i,1, 1.0 );
    }
}
// ***********************************************
