#include <cmath>
#include "DogParams.h"
#include "tensors.h"
#include "dog_math.h"
#include "stdio.h"
#include <cstdlib>
using namespace std;

// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x = Psi(q,x,t)
//
void ConstructL_Interface_ImplicitPart(const int method[],
        const dTensor2& node,
        dTensorBC3& aux,
        dTensorBC3& q,      // setbndy conditions modifies q
        dTensor4 qI,
        dTensor4 auxI,
        dTensorBC1 global2interf,
        dTensor1 interf2global,
        dTensor2 dxi,
        dTensor4& LstarI,
        dTensorBC3& Lstar,dTensor4& Implicit,
        dTensorBC1& smax)
{

void L2Project_interface(int mopt, int istart, int iend,
               const dTensor2& node,
               const dTensorBC3& qin,
               const dTensorBC3& auxin,
               const dTensor4 qI,
               const dTensor4 auxI,dTensor4& Implicit,double dt,
               const dTensor1 interf2global,
               const dTensorBC1 global2interf,
               const dTensor2 dxi,
               dTensorBC3& Fout,
               dTensor4& FI,
               void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&));

    double RiemannSolve(const dTensor1& xedge,
            const dTensor1& Ql,
            const dTensor1& Qr,
            const dTensor1& Auxl,
            const dTensor1& Auxr,
            dTensor1& Fl,
            dTensor1& Fr,
            void (*FluxFunc)(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&),
            void (*SetWaveSpd)(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
                const dTensor1&,double&,double&));
    void SetBndValues(const dTensor2&,dTensorBC3&,dTensorBC3&);
    void L2Project(int,int,int,const dTensor2&,const dTensorBC3&,const dTensorBC3&,dTensorBC3&,
            void (*Func)(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&));
    void FluxFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void SetWaveSpd(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
            const dTensor1&,double&,double&);
    void SourceTermFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void LstarExtra(const dTensor2&,dTensorBC3&,dTensorBC3&,dTensorBC3&);
    void ArtificialViscosity(const dTensor2&,dTensorBC3&,dTensorBC3&,dTensorBC3&);

    const int melems = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   kmax = q.getsize(3);
    const int   maux = aux.getsize(2);
    const int    mbc = q.getmbc();
    const int   interfaces = interf2global.getsize();

     
    // Flux values, interior integrals and sourse term, respectively.
    dTensorBC2 Fm(melems,meqn,mbc);
    dTensorBC2 Fp(melems,meqn,mbc);
    dTensorBC3  N(melems,meqn,kmax,mbc);
    dTensor4  NI(2,interfaces,meqn,kmax);NI.setall(0.0);
    dTensorBC3 Psi(melems,meqn,kmax,mbc);

    dTensor3  Finterface(interfaces,meqn,2);Finterface.setall(0.0);



    // Grid spacing -- node( 1:(mx+1), 1 ) = cell edges
    const double xlower = node.get(1,1);
    const double     dx = node.get(2,1) - node.get(1,1);

    // ---------------------------------------------------------
    // Part I: compute inter-element interaction fluxes
    // ---------------------------------------------------------

    // Boundary conditions
    SetBndValues(node,aux,q);

    // Loop over interior edges and solve Riemann problems
    //#pragma omp parallel for
    //
    // This sets both smax(i) and smax(i-1), so can't be parallelized!
    //

    int jlist[4];
    jlist[0]=-2;
    jlist[1]=-1;
    jlist[2]=1;
    jlist[3]=1;

    int startlist[4];
    startlist[0]=0;
    startlist[1]=kmax;
    startlist[2]=4*kmax;
    startlist[3]=5*kmax;

    for (int i=(2-mbc); i<=(melems+mbc); i++)
    {

      if(abs(global2interf.get(i))<=1.0e-1)
      {
      }
      else if(abs(global2interf.get(i))>1.0e-1)
      {

       int iint=int(global2interf.get(i));

       for(int j1=0;j1<4;j1++)
       {
        int j=jlist[j1];

        dTensor1   Ql(meqn),   Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);

        // Riemann data - q
        for (int m=1; m<=meqn; m++)
        {

            for (int k=1; k<=method[1]; k++)
            for (int k1=1; k1<=method[1]; k1++)
            {
                double tmpl=-pow(-1.0,k-1)*pow(-1.0,k1-1)*sqrt(2.0*double(k)-1.0)*sqrt(2.0*double(k1)-1.0)*aux.get(i+j,1,1);
                double tmpr=pow(1.0,k-1)*pow(1.0,k1-1)*sqrt(2.0*double(k)-1.0)*sqrt(2.0*double(k1)-1.0)*aux.get(i+j,1,1);
                if(j==-2)
                {Implicit.set(iint,m,startlist[j+2]+k,startlist[j+2]+k1,tmpl+tmpr);}
                else
                {
                    Implicit.set(iint,m,startlist[j+2]+k,startlist[j+1]+k1,tmpl);
                    Implicit.set(iint,m,startlist[j+2]+k,startlist[j+2]+k1,tmpr);
                }
            }
        }

      
       }
     
       
 
       double xcs[2];
       double xc1 = xlower + (double(i)-1.0)*dx+0.5*dxi.get(iint,1);
       double xc2 = xlower + (double(i)-1.0)*dx+dxi.get(iint,1)+0.5*dxi.get(iint,2);
 
       double xc = xlower + (double(i) - 0.5)*dx;
       xcs[0]=xc1;
       xcs[1]=xc2;
       
       dTensor1 Qlocal(kmax);
       int mtmp = iMax(1,6);//iMax(1,mpoints-mopt);
       dTensor1 xpts(mtmp);
       dTensor2 qvalsl(mtmp,meqn);
       dTensor2 qvalsr(mtmp,meqn);
       dTensor1 wgt(mtmp), spts(mtmp);

      dTensor3 phiI(mtmp,kmax,2);
       dTensor2 phi(mtmp,kmax);
      // Set quadrature weights and points
      void SetQuadPts(dTensor1&,dTensor1&);
      SetQuadPts(wgt,spts);
      // Sample basis function at quadrature points
      void SampleBasis_Interface(const dTensor1&,dTensor3&,double xc,double dx,double x1,double x2,double dx1,double dx2);
      SampleBasis_Interface(spts,phiI,xc,dx,xc1,xc2,dxi.get(iint,1),dxi.get(iint,2));

      void SampleBasis(const dTensor1&,dTensor2&);
      SampleBasis(spts,phi);

      for(int i1=1;i1<=2;i1++)
       { dTensor1   Ql(meqn),   Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);

        // Riemann data - q
        for (int m=1; m<=meqn; m++)
        {

            if(i1==1)
            {for (int k=1; k<=kmax;k++)//method[1]; k++)
            for (int k1=1; k1<=kmax;k1++)
             {

                double tmpl=-pow(-1.0,k-1)*pow(-1.0,k1-1)*sqrt(2.0*double(k)-1.0)*sqrt(2.0*double(k1)-1.0)*auxI.get(i1,iint,1,1);
                double tmpr=pow(1.0,k-1)*pow(1.0,k1-1)*sqrt(2.0*double(k)-1.0)*sqrt(2.0*double(k1)-1.0)*auxI.get(i1,iint,1,1);

                Implicit.set(iint,m,2*kmax+k,kmax+k1,tmpl);
                Implicit.set(iint,m,2*kmax+k,2*kmax+k1,tmpr);
            }
            }

            if(i1==2)
            {
            for (int k=1; k<=kmax;k++)//method[1]; k++)
            for (int k1=1; k1<=kmax;k1++)
            {
                double tmpl=-pow(-1.0,k-1)*pow(-1.0,k1-1)*sqrt(2.0*double(k)-1.0)*sqrt(2.0*double(k1)-1.0)*auxI.get(i1,iint,1,1);
                double tmpr=pow(1.0,k-1)*pow(1.0,k1-1)*sqrt(2.0*double(k)-1.0)*sqrt(2.0*double(k1)-1.0)*auxI.get(i1,iint,1,1);

                    Implicit.set(iint,m,3*kmax+k,2*kmax+k1,tmpl);
                    Implicit.set(iint,m,3*kmax+k,3*kmax+k1,tmpr);
            }
            }
        }
       
        // Solve Riemann problem
        dTensor1 xedge(1);
        if(i1==1)
        {xedge.set(1, xlower + (double(i)-1.0)*dx );}
        if(i1==2)
        {xedge.set(1,xlower + (double(i)-1.0)*dx+dxi.get(iint,1));}

        dTensor1 Fl(meqn);
        dTensor1 Fr(meqn);
        double smax_edge = RiemannSolve(xedge,Ql,Qr,Auxl,Auxr,Fl,Fr,
                &FluxFunc,&SetWaveSpd);

        // This is a problem for the pragma statements! (-DS)
        smax.set(i-1, Max(smax_edge,smax.get(i-1)) );
        smax.set(i,   Max(smax_edge,smax.get(i)) );
        

      }
  }
  }
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part II: compute intra-element contributions
    // ---------------------------------------------------------
    //
    //   N = int( f(q,x,t) * phi_x, x )/dx
    //
    // Compute ``N'' by projecting flux function onto the 
    // gradient of Legendre polynomials
    if (method[1]>1)
    {  
       //L2Project(1,1-mbc,melems+mbc,node,q,aux,N,&FluxFunc);
       L2Project_interface(1,1-mbc,melems+mbc,node,q,aux,qI,auxI,Implicit,0.0,interf2global,global2interf,dxi,N,NI,&FluxFunc);
    }
    else
    {
#pragma omp parallel for
        for (int i=1-mbc; i<=(melems+mbc); i++)
            for (int m=1; m<=meqn; m++)
            {
                N.set(i,m,1, 0.0 );   
            }	
    }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part III: compute source term
    // --------------------------------------------------------- 
    if ( method[7]>0 )
    {        
        // Set source term on computational grid
        // Set values and apply L2-projection
        L2Project(0,1-mbc,melems+mbc,node,q,aux,Psi,&SourceTermFunc);
    }
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part IV: construct Lstar
    // ---------------------------------------------------------
    /*if (method[7]==0)  // Without Source Term
    {
#pragma omp parallel for
        for (int i=(2-mbc); i<=(melems+mbc-1); i++)	
            for (int m=1; m<=meqn; m++)
                for (int k=1; k<=method[1]; k++)
                {
                    if(abs(global2interf.get(i))<=1.0e-1)
                    {
                    double tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*
                        ( Fp.get(i,m) + pow(-1.0,k)*Fm.get(i,m) )/dx;

                    Lstar.set(i,m,k, tmp );
                    }
                    else
                    {
                        
                        int iint=int(global2interf.get(i)); 

                        //double tmp = NI.get(1,iint,m,k)-sqrt(2.0*double(k)-1.0)*(Finterface.get(iint,m,1)+pow(-1.0,k)*Fm.get(i,m))/dxi.get(iint,1);
                        double tmp = -sqrt(2.0*double(k)-1.0)*(Finterface.get(iint,m,1)+pow(-1.0,k)*Fm.get(i,m));
                        LstarI.set(1,iint,m,k,tmp);
                            //if(abs(LstarI.get(1,iint,m,k))>1.0e-4)
                         
                        //tmp = NI.get(2,iint,m,k)-sqrt(2.0*double(k)-1.0)*(Fp.get(i,m)+pow(-1.0,k)*Finterface.get(iint,m,2))/dxi.get(iint,2);
                        tmp = -sqrt(2.0*double(k)-1.0)*(Fp.get(i,m)+pow(-1.0,k)*Finterface.get(iint,m,2));
                        LstarI.set(2,iint,m,k,tmp);
                    }	      
                }
    }
    else  // With Source Term
    {
#pragma omp parallel for
        for (int i=(2-mbc); i<=(melems+mbc-1); i++)	
            for (int m=1; m<=meqn; m++)	    
                for (int k=1; k<=method[1]; k++)
                {
                    double tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*
                        ( Fp.get(i,m) + pow(-1.0,k)*Fm.get(i,m) )/dx
                        + Psi.get(i,m,k);

                    Lstar.set(i,m,k, tmp );
                }
    }*/
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part V: add extra contributions to Lstar
    // ---------------------------------------------------------
    // Call LstarExtra
    LstarExtra(node,aux,q,Lstar);
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part VI: artificial viscosity limiter
    // ---------------------------------------------------------
    if (dogParams.get_space_order()>1)
    {
        if (dogParams.using_viscosity_limiter())
        {  ArtificialViscosity(node,aux,q,Lstar);  }
    }
    // ---------------------------------------------------------

}


