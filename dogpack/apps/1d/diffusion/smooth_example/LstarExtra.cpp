#include "tensors.h"
#include "diffusion.h"

// Compute explicit time-stepping diffusion operator
void LstarExtra(const dTensor2& node, dTensorBC3& aux, 
        dTensorBC3& q, dTensorBC3& Lstar)
{
    int i,k,ell;
    double x; 
    double Stmp;
    const int    mx = q.getsize(1);
    const int  meqn = q.getsize(2);
    const int  kmax = q.getsize(3);
    const int   mbc = q.getmbc();
    const int  maux = aux.getsize(2);
    const double dx = node.get(2,1)-node.get(1,1);
    dTensorBC3  V(mx,meqn,kmax,mbc);
    dTensorBC2 FU(mx,1,mbc);
    dTensorBC2 FV(mx,1,mbc);

    // Set boundary conditions
    void SetBndValues(const dTensor2&,dTensorBC3&,dTensorBC3&);
    SetBndValues(node,aux,q);

    // fluxes to evaluate gradient
    for (i=(1-mbc); i<=(mx+mbc); i++)
    {
        FU.set(i,1, -signv[0]*phiv[0]*q.get(i,1,1) );

        for (k=2; k<=kmax; k++)
        {  FU.set(i,1, FU.get(i,1) - signv[k-1]*phiv[k-1]*q.get(i,1,k) );  }
    }

    // gradient
    for (i=(1-mbc); i<=(mx+mbc-1); i++)
        for (k=1; k<=kmax; k++)
        {
            Stmp = Smat[k-1][0]*q.get(i,1,1);
            for (ell=2; ell<=kmax; ell++)
            {  Stmp = Stmp + Smat[k-1][ell-1]*q.get(i,1,ell);  }

            V.set(i,1,k, (phiv[k-1]*(FU.get(i+1,1) + 
                            signv[k-1]*FU.get(i,1)) - Stmp)/dx );
        }

    // fluxes to evaluate solution
    for (i=(2-mbc); i<=(mx+mbc); i++)
    {
        FV.set(i,1, phiv[0]*V.get(i-1,1,1) );

        for (k=2; k<=kmax; k++)
        {
            FV.set(i,1, FV.get(i,1) + phiv[k-1]*V.get(i-1,1,k) );      
        }
    }

    // update solution
    for (i=(2-mbc); i<=(mx+mbc-1); i++)
        for (k=1; k<=kmax; k++)
        {
            Stmp = Smat[k-1][0]*V.get(i,1,1);
            for (ell=2; ell<=kmax; ell++)
            {  Stmp = Stmp + Smat[k-1][ell-1]*V.get(i,1,ell);  }

            Lstar.set(i,1,k, Lstar.get(i,1,k) 
                    + (phiv[k-1]*(FV.get(i+1,1) + 
                            signv[k-1]*FV.get(i,1)) - Stmp)/dx );
        }
}
