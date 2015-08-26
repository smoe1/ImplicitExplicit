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
void ConstructL_Interface_ImplicitBox(const int method[],
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
    for (int i=(2-mbc); i<=(melems+mbc); i++)
    {

      if(abs(global2interf.get(i))<=1.0e-1 && (abs(global2interf.get(i-1))<=1.0e-1))
      {
        dTensor1   Ql(meqn),   Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);

        // Riemann data - q
        for (int m=1; m<=meqn; m++)
        {
            Ql.set(m, 0.0 );
            Qr.set(m, 0.0 );

            for (int k=1; k<=method[1]; k++)
            {
                Ql.set(m, Ql.get(m) + sqrt(2.0*double(k)-1.0)
                        *q.get(i-1,m,k) );
                Qr.set(m, Qr.get(m) + pow(-1.0,k+1)
                        *sqrt(2.0*double(k)-1.0)*q.get(i,m,k) ); 
            }
        }

        // Riemann data - aux
        for (int m=1; m<=maux; m++)
        {
            Auxl.set(m, 0.0 );
            Auxr.set(m, 0.0 );

            for (int k=1; k<=method[1]; k++)
            {
                Auxl.set(m, Auxl.get(m) + sqrt(2.0*double(k)-1.0)
                        *aux.get(i-1,m,k) );
                Auxr.set(m, Auxr.get(m) + pow(-1.0,k+1)
                        *sqrt(2.0*double(k)-1.0)*aux.get(i,m,k) ); 
            }
        }

        // Solve Riemann problem
        dTensor1 xedge(1);
        xedge.set(1, xlower + (double(i)-1.0)*dx );

        dTensor1 Fl(meqn);
        dTensor1 Fr(meqn);
        double smax_edge = RiemannSolve(xedge,Ql,Qr,Auxl,Auxr,Fl,Fr,
                &FluxFunc,&SetWaveSpd);

        // This is a problem for the pragma statements! (-DS)
        smax.set(i-1, Max(smax_edge,smax.get(i-1)) );
        smax.set(i,   Max(smax_edge,smax.get(i)) );

        // Construct fluxes
        for (int m=1; m<=meqn; m++)
        {
            Fm.set(i,  m, Fr.get(m) );
            Fp.set(i-1,m, Fl.get(m) );
        }

      }
      else if(abs(global2interf.get(i))>1.0e-1)
      {

       int iint=int(global2interf.get(i));
 
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

          int mpoints=kmax;
          for (int m=1; m<=mtmp; m++)
            {
              // Solution values (q) at each grid point
              for (int me=1; me<=meqn; me++)
                {
                  qvalsl.set(m,me, 0.0 );
                  qvalsr.set(m,me, 0.0 );
                  for (int k=1; k<=mpoints; k++)
                    {
                      qvalsl.set(m,me, qvalsl.get(m,me)
                                + phi.get(m,k) * qI.get(1,iint,me,k) );
                      qvalsr.set(m,me, qvalsr.get(m,me)
                                + phi.get(m,k) * qI.get(2,iint,me,k) );
                    }
                }
            }
       for (int k1=1;k1<=kmax;k1++)
       {
          Qlocal.set(k1,0.0);
          double Qi=0.0;
          double tmpl=0.0;
          double tmpr=0.0;
          for (int k=1; k<=mtmp; k++)
          {
              tmpl += wgt.get(k)*qvalsl.get(k,1)*phiI.get(k,k1,1);
              tmpr += wgt.get(k)*qvalsr.get(k,1)*phiI.get(k,k1,2);
          }
          Qlocal.set(k1,(0.5*dxi.get(iint,1)*tmpl+0.5*dxi.get(iint,2)*tmpr)/dx);
       }

       for(int i1=1;i1<=2;i1++)
       { dTensor1   Ql(meqn),   Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);

        // Riemann data - q
        for (int m=1; m<=meqn; m++)
        {

            Ql.set(m, 0.0 );
            Qr.set(m, 0.0 );
            if(i1==1)
            {for (int k=1; k<=kmax;k++)//method[1]; k++)
            {
                double average=(dxi.get(iint,i1)*qI.get(i1,iint,m,k)+dxi.get(iint,i1+1)*qI.get(i1+1,iint,m,k))/dx;
                Ql.set(m, Ql.get(m) + sqrt(2.0*double(k)-1.0)
                        *q.get(i-1,m,k) );
                //Qr.set(m, Qr.get(m) + pow(-1.0,k+1)
                //        *sqrt(2.0*double(k)-1.0)*qI.get(i1,iint,m,k) ); 
                Qr.set(m, Qr.get(m) + pow(-1.0,k+1)
                        *sqrt(2.0*double(k)-1.0)*qI.get(i1,iint,m,k)  ); 
                printf("%d average value = %e %e %e %e \n",k,average,qI.get(i1,iint,m,k),qI.get(i1+1,iint,m,k),dxi.get(iint,i1)+dxi.get(iint,i1+1)-dx);
                //Qr.set(m, Qr.get(m) + pow(-1.0,k+1)*sqrt(2.0*double(k)-1.0)*average ); 
            }
            }

            if(i1==2)
            {for (int k=1; k<=kmax;k++)//method[1]; k++)
            {
                Ql.set(m, Ql.get(m) + sqrt(2.0*double(k)-1.0)
                        *qI.get(i1-1,iint,m,k) );
                Qr.set(m, Qr.get(m) + pow(-1.0,k+1)
                        *sqrt(2.0*double(k)-1.0)*qI.get(i1,iint,m,k) ); 
            }
            }
        }
       
        // Riemann data - aux
        for (int m=1; m<=maux; m++)
        {
            Auxl.set(m, 0.0 );
            Auxr.set(m, 0.0 );
            if(i1==1)
            {
            for (int k=1; k<=method[1]; k++)
            {
                Auxl.set(m, Auxl.get(m) + sqrt(2.0*double(k)-1.0)
                        *aux.get(i-1,m,k) );
                Auxr.set(m, Auxr.get(m) + pow(-1.0,k+1)
                        *sqrt(2.0*double(k)-1.0)*auxI.get(i1,iint,m,k) ); 
            }
            //printf("Aux= %e %e \n",Auxl.get(1),Auxr.get(1));
            }
            if(i1==2)
            {
            for (int k=1; k<=method[1]; k++)
            {
                Auxl.set(m, Auxl.get(m) + sqrt(2.0*double(k)-1.0)
                        *auxI.get(i1-1,iint,m,k) );
                Auxr.set(m, Auxr.get(m) + pow(-1.0,k+1)
                        *sqrt(2.0*double(k)-1.0)*auxI.get(i1,iint,m,k) );
            }
            }
        }
        // Solve Riemann problem
        dTensor1 xedge(1);
        if(i1==1)
        {xedge.set(1, xlower + (double(i)-1.0)*dx );}
        if(i1==2)
        {xedge.set(1,xlower + (double(i)-1.0)*dx+dxi.get(iint,1));

         //printf("x=%e Aux= %e %e \n",xedge.get(1),Auxl.get(1),Auxr.get(1));
        }

        dTensor1 Fl(meqn);
        dTensor1 Fr(meqn);
        double smax_edge = RiemannSolve(xedge,Ql,Qr,Auxl,Auxr,Fl,Fr,
                &FluxFunc,&SetWaveSpd);

        // This is a problem for the pragma statements! (-DS)
        smax.set(i-1, Max(smax_edge,smax.get(i-1)) );
        smax.set(i,   Max(smax_edge,smax.get(i)) );
        


        // Construct fluxes
        for (int m=1; m<=meqn; m++)
        {
            if(i1==2)
            {
             //Finterface.set(iint,m,1,Fl.get(m));
             //Finterface.set(iint,m,2,Fr.get(m));
             Finterface.set(iint,m,1,0.0);
             Finterface.set(iint,m,2,0.0);

             //printf("HERE %e %e %e %e \n",Fl.get(m),Fr.get(m),qI.get(1,1,1,1),qI.get(2,1,1,1));
            }

            if(i1==1)
            {
            //Fm.set(i,  m, Fr.get(m) );
            //Fp.set(i-1,m, Fl.get(m) );
            Fm.set(i,  m, 0.0 );
            Fp.set(i-1,m, 0.0 );
            
            double Us1=auxI.get(1,iint,1,1);
            double Us2=auxI.get(2,iint,1,1);
       
            printf("here %e %e \n",Fm.get(i,1),Fp.get(i-1,1));

            
 
	    int k=1;
             for(int il=1;il<=5*kmax;il++)
             for(int jl=1;jl<=5*kmax;jl++)
             {Implicit.set(iint,1,il,jl,0.0);}
            if(kmax==1)
            {
               Implicit.set(iint,1,1,1,Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(1)-1.0));             
               Implicit.set(iint,1,2,1,-Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(1)-1.0));             
               Implicit.set(iint,1,2,2,Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(1)-1.0));             
               Implicit.set(iint,1,3,2,-Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(1)-1.0));             
               Implicit.set(iint,1,3,3,Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(1)-1.0));             
               Implicit.set(iint,1,4,3,Us2*sqrt(2.0*double(1)-1.0)*pow(-1.0,1)*sqrt(2.0*double(1)-1.0));
               Implicit.set(iint,1,4,4,Us2*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(1)-1.0));
               Implicit.set(iint,1,5,4,Us2*sqrt(2.0*double(1)-1.0)*pow(-1.0,1)*sqrt(2.0*double(1)-1.0));
               Implicit.set(iint,1,5,5,Us2*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(1)-1.0));
               //Implicit.set(iint,1,6,5,Us2*sqrt(2.0*double(1)-1.0)*pow(-1.0,1)*sqrt(2.0*double(1)-1.0));
            }
            if(kmax==2)
            {

            Implicit.set(iint,1,1,1,Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,1,2,Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(2)-1.0));             
            Implicit.set(iint,1,2,1,Us1*sqrt(2.0*double(2)-1.0)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,2,2,Us1*sqrt(2.0*double(2)-1.0)*sqrt(2.0*double(2)-1.0));

            Implicit.set(iint,1,3,1,Us2*sqrt(2.0*double(1)-1.0)*pow(-1.0,1)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,3,2,Us2*sqrt(2.0*double(1)-1.0)*pow(-1.0,1)*sqrt(2.0*double(2)-1.0));
            Implicit.set(iint,1,4,1,Us2*sqrt(2.0*double(2)-1.0)*pow(-1.0,2)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,4,2,Us2*sqrt(2.0*double(2)-1.0)*pow(-1.0,2)*sqrt(2.0*double(2)-1.0));

            Implicit.set(iint,1,3,3,Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,3,4,Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(2)-1.0));
            Implicit.set(iint,1,4,3,Us1*sqrt(2.0*double(2)-1.0)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,4,4,Us1*sqrt(2.0*double(2)-1.0)*sqrt(2.0*double(2)-1.0));

            Implicit.set(iint,1,5,3,Us2*sqrt(2.0*double(1)-1.0)*pow(-1.0,1)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,5,4,Us2*sqrt(2.0*double(1)-1.0)*pow(-1.0,1)*sqrt(2.0*double(2)-1.0));
            Implicit.set(iint,1,6,3,Us2*sqrt(2.0*double(2)-1.0)*pow(-1.0,2)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,6,4,Us2*sqrt(2.0*double(2)-1.0)*pow(-1.0,2)*sqrt(2.0*double(2)-1.0));

            Implicit.set(iint,1,5,5,Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,5,6,Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(2)-1.0));
            Implicit.set(iint,1,6,5,Us1*sqrt(2.0*double(2)-1.0)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,6,6,Us1*sqrt(2.0*double(2)-1.0)*sqrt(2.0*double(2)-1.0));

            Implicit.set(iint,1,7,5,Us2*sqrt(2.0*double(1)-1.0)*pow(-1.0,1)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,7,6,Us2*sqrt(2.0*double(1)-1.0)*pow(-1.0,1)*sqrt(2.0*double(2)-1.0));
            Implicit.set(iint,1,8,5,Us2*sqrt(2.0*double(2)-1.0)*pow(-1.0,2)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,8,6,Us2*sqrt(2.0*double(2)-1.0)*pow(-1.0,2)*sqrt(2.0*double(2)-1.0));

            Implicit.set(iint,1,7,7,Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,7,8,Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(2)-1.0));
            Implicit.set(iint,1,8,7,Us1*sqrt(2.0*double(2)-1.0)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,8,8,Us1*sqrt(2.0*double(2)-1.0)*sqrt(2.0*double(2)-1.0));

            Implicit.set(iint,1,9,7,Us2*sqrt(2.0*double(1)-1.0)*pow(-1.0,1)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,9,8,Us2*sqrt(2.0*double(1)-1.0)*pow(-1.0,1)*sqrt(2.0*double(2)-1.0));
            Implicit.set(iint,1,10,7,Us2*sqrt(2.0*double(2)-1.0)*pow(-1.0,2)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,10,8,Us2*sqrt(2.0*double(2)-1.0)*pow(-1.0,2)*sqrt(2.0*double(2)-1.0));

            Implicit.set(iint,1,9,9,Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,9,10,Us1*sqrt(2.0*double(1)-1.0)*sqrt(2.0*double(2)-1.0));
            Implicit.set(iint,1,10,9,Us1*sqrt(2.0*double(2)-1.0)*sqrt(2.0*double(1)-1.0));
            Implicit.set(iint,1,10,10,Us1*sqrt(2.0*double(2)-1.0)*sqrt(2.0*double(2)-1.0));




            }            
            }
        }
      }  
      }
      else if(abs(global2interf.get(i-1))>1.0e-1)
      {
       int iint=int(global2interf.get(i-1));
       int i1=2;
        dTensor1   Ql(meqn),   Qr(meqn);
        dTensor1 Auxl(maux), Auxr(maux);


       double xc1 = xlower + (double(i-1)-1.0)*dx+0.5*dxi.get(iint,1);
       double xc2 = xlower + (double(i-1)-1.0)*dx+dxi.get(iint,1)+0.5*dxi.get(iint,2);
       double xc = xlower + (double(i-1) - 0.5)*dx;
       
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

      int mpoints=kmax;
          for (int m=1; m<=mtmp; m++)
            {
              // Solution values (q) at each grid point
              for (int me=1; me<=meqn; me++)
                {
                  qvalsl.set(m,me, 0.0 );
                  qvalsr.set(m,me, 0.0 );
                  for (int k=1; k<=mpoints; k++)
                    {
                      qvalsl.set(m,me, qvalsl.get(m,me)
                                + phi.get(m,k) * qI.get(1,iint,me,k) );
                      qvalsr.set(m,me, qvalsr.get(m,me)
                                + phi.get(m,k) * qI.get(2,iint,me,k) );
                    }
                }
            }
       for (int k1=1;k1<=kmax;k1++)
       {
          Qlocal.set(k1,0.0);
          double Qi=0.0;
          double tmpl=0.0;
          double tmpr=0.0;
          for (int k=1; k<=mtmp; k++)
          {
              tmpl += wgt.get(k)*qvalsl.get(k,1)*phiI.get(k,k1,1);
              tmpr += wgt.get(k)*qvalsr.get(k,1)*phiI.get(k,k1,2);
          }
          Qlocal.set(k1,0.5*(dxi.get(iint,1)*tmpl+dxi.get(iint,2)*tmpr)/dx);
       }
       
        //printf("%d HERE! %e %e\n",i1,xc1+dxi.get(iint,1)/2.0,xc-dx/2.0);

        // Riemann data - q
        for (int m=1; m<=meqn; m++)
        {
            Ql.set(m, 0.0 );
            Qr.set(m, 0.0 );

            
            {for (int k=1; k<=kmax; k++)
            {
                double average=(dxi.get(iint,1)*qI.get(1,iint,m,k)+dxi.get(iint,2)*qI.get(2,iint,m,k))/(dx);
                //printf("%d average value2 = %e %e %e %e \n",k,average,qI.get(i1-1,iint,m,k),qI.get(i1,iint,m,k),dxi.get(iint,i1-1)+dxi.get(iint,i1)-dx);

                Ql.set(m, Ql.get(m) + sqrt(2.0*double(k)-1.0)
                        *Qlocal.get(k) );
                //Ql.set(m, Ql.get(m) + sqrt(2.0*double(k)-1.0)*average );
                //Ql.set(m, Ql.get(m) + sqrt(2.0*double(k)-1.0)
                //        *qI.get(2,iint,m,k) );
                Qr.set(m, Qr.get(m) + pow(-1.0,k+1)
                        *sqrt(2.0*double(k)-1.0)*q.get(i,m,k) ); 
            }
            }
        }
        // Riemann data - aux
        for (int m=1; m<=maux; m++)
        {
            Auxl.set(m, 0.0 );
            Auxr.set(m, 0.0 );
            for (int k=1; k<=method[1]; k++)
            {
                Auxl.set(m, Auxl.get(m) + sqrt(2.0*double(k)-1.0)*auxI.get(i1,iint,m,k) );
                Auxr.set(m, Auxr.get(m) + pow(-1.0,k+1)
                        *sqrt(2.0*double(k)-1.0)*aux.get(i,m,k) );
            }            
        }
        // Solve Riemann problem
        dTensor1 xedge(1);
        xedge.set(1, xlower + (double(i)-1.0)*dx );
        dTensor1 Fl(meqn);
        dTensor1 Fr(meqn);
        double smax_edge = RiemannSolve(xedge,Ql,Qr,Auxl,Auxr,Fl,Fr,
                &FluxFunc,&SetWaveSpd);
        // This is a problem for the pragma statements! (-DS)
        smax.set(i-1, Max(smax_edge,smax.get(i-1)) );
        smax.set(i,   Max(smax_edge,smax.get(i)) );
        // Construct fluxes
           for (int m=1; m<=meqn; m++)
           {
              //Fm.set(i,  m, Fr.get(m) );
              //Fp.set(i-1,m, Fl.get(m) );
              Fm.set(i,  m, 0.0 );
              Fp.set(i-1,m, 0.0 );
           }
      }
      else{printf("Warning! Should not be here! \n");exit(1);}
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
    if (method[7]==0)  // Without Source Term
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
    }
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


