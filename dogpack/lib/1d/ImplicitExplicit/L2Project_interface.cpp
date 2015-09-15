#include "dogdefs.h"
#include "stdlib.h"
#include "DogParams.h"
#include "DogParamsCart1.h"
#include "dog_math.h"

// All-purpose routine for computing the L2-projection
// of various functions onto:
//     mopt==0:   the Legendre basis
//     mopt==1:   the derivatives of Legendre basis
//
// ---------------------------------------------------------------------
// Inputs should have the following sizes:   
//           dTensor2    node(mnodes,1)
//           dTensorBC3 auxin(1-mbc:mnodes+mbc,maux,mpoints)
//           dTensorBC3   qin(1-mbc:mnodes+mbc,meqn,mpoints)
//           dTensorBC3  Fout(1-mbc:mnodes+mbc,mlength,mpoints)
// ---------------------------------------------------------------------

void L2Project_interface(int mopt, int run,int istart, int iend,
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
	       void (*Func)(const dTensor1&, const dTensor2&, const dTensor2&, dTensor2&))
{    
  const int kmax = dogParams.get_space_order();
  const int meqn = qin.getsize(2);
  const int maux = auxin.getsize(2);
  const int mlength = Fout.getsize(2);
  const int mpoints = Fout.getsize(3);  

  int mtmp = iMax(1,6-mopt);//iMax(1,mpoints-mopt);
  dTensor1 wgt(mtmp), spts(mtmp),rspts(mtmp),lspts(mtmp);
  dTensor2 phi(mtmp,kmax), phi_x(mtmp,kmax), rIphi_x(mtmp,kmax), lIphi_x(mtmp,kmax);

  // -----------------
  // Quick error check
  // -----------------
  if (meqn<1 || maux <0 || mpoints<1 || mpoints>6 || mlength<1 
      || mopt < 0 || mopt > 1)
    {
      printf(" Error in L2project.cpp ... \n");
      printf("         meqn = %i\n",meqn);
      printf("         maux = %i\n",maux);
      printf("      mpoints = %i\n",mpoints);
      printf("      mlength = %i\n",mlength);
      printf("       istart = %i\n",istart);
      printf("         iend = %i\n",iend);
      printf("        mopts = %i\n",mopt);
      printf("\n");
      exit(1);
    }

  // ---------------------------------------------
  // Check for trivial case in the case of mopt==1
  // ---------------------------------------------
  if ( mpoints == mopt )
    { Fout.setall(0.); }
  else
    {
      // Set quadrature weights and points
      void SetQuadPts(dTensor1&,dTensor1&);
      SetQuadPts(wgt,spts);

      // Sample basis function at quadrature points
      void SampleBasis(const dTensor1&,dTensor2&);
      SampleBasis(spts,phi);

      // Sample gradient of basis function at quadrature points
      void SampleBasisGrad(const dTensor1&,dTensor2&);
      void SampleBasisGrad_variable(const dTensor1&,dTensor2&,double);
      SampleBasisGrad(spts,phi_x);

      // ----------------------------------
      // Loop over all elements of interest
      // ----------------------------------    
      const double xlow = dogParamsCart1.get_xlow();
      const double dx = dogParamsCart1.get_dx();

#pragma omp parallel for
      for (int i=istart; i<=iend; i++)
        {
          if(abs(global2interf.get(i))<1.0e-1)
          {
	  double xc = xlow + (double(i)-0.5)*dx;

	  // Each of these three items needs to be private to each thread ..
	  dTensor1 xpts(mtmp);
	  dTensor2 qvals(mtmp,meqn);
	  dTensor2 auxvals(mtmp,maux);
	  dTensor2 fvals(mtmp,mlength);

	  // Loop over each quadrature point
	  for (int m=1; m<=mtmp; m++)
            {
	      // grid point x
	      xpts.set( m, xc + 0.5*dx*spts.get(m) );

	      // Solution values (q) at each grid point
	      for (int me=1; me<=meqn; me++)
                {
		  qvals.set(m,me, 0.0 );

		  for (int k=1; k<=mpoints; k++)
                    {
		      qvals.set(m,me, qvals.get(m,me) 
                                + phi.get(m,k) * qin.get(i,me,k) );
                    }
                }

	      // Auxiliary values (aux) at each grid point
	      for (int ma=1; ma<=maux; ma++)
                {
		  auxvals.set(m,ma, 0.0 );

		  for (int k=1; k<=mpoints; k++)
                    {
		      auxvals.set(m,ma, auxvals.get(m,ma) 
				  + phi.get(m,k) * auxin.get(i,ma,k) );
                    }
                }
            }


	  // Call user-supplied function to set fvals
	  Func(xpts,qvals,auxvals,fvals);

	  // Evaluate integrals
	  if (mopt==0) // project onto Legendre basis
            {
	      for (int m1=1; m1<=mlength; m1++)
		for (int m2=1; m2<=mpoints; m2++)
		  {
		    double tmp = 0.0;
		    for (int k=1; k<=mtmp; k++)
		      {
			tmp += wgt.get(k)*fvals.get(k,m1)*phi.get(k,m2);
		      }
		    Fout.set(i,m1,m2, 0.5*tmp );
		  }

            }
	  else // project onto derivatives of Legendre basis
            {
	      for (int m1=1; m1<=mlength; m1++)             
		for (int m2=1; m2<=mpoints; m2++)
		  {
		    double tmp = 0.0;
		    for (int k=1; k<=mtmp; k++)
		      {
			tmp += wgt.get(k)*fvals.get(k,m1)*phi_x.get(k,m2);
		      }
		    Fout.set(i,m1,m2, 0.5*tmp );
		  }
            }
        }
          else if(run==1)
          {


          int iint=int(global2interf.get(i));

	  double xc1 = xlow + (double(i)-1.0)*dx+0.5*dxi.get(iint,1);
          double xc2 = xlow + (double(i)-1.0)*dx+dxi.get(iint,1)+0.5*dxi.get(iint,2);

          double dxl=dxi.get(iint,1);
          double dxr=dxi.get(iint,2);

        
	  // Each of these three items needs to be private to each thread ..
	  dTensor1 xptsl(mtmp);
	  dTensor2 qvalsl(mtmp,meqn);
	  dTensor2 auxvalsl(mtmp,maux);
	  dTensor2 fvalsl(mtmp,mlength);


	  dTensor1 xptsr(mtmp);
	  dTensor2 qvalsr(mtmp,meqn);
	  dTensor2 auxvalsr(mtmp,maux);
	  dTensor2 fvalsr(mtmp,mlength);

          //SampleBasisGrad_variable(lspts,lIphi_x,dxl);
          //SampleBasisGrad_variable(rspts,rIphi_x,dxr);           
          //SampleBasisGrad_variable(lspts,lIphi_x,1.0);
          //SampleBasisGrad_variable(rspts,rIphi_x,1.0);           

          SampleBasisGrad_variable(spts,lIphi_x,1.0);
          SampleBasisGrad_variable(spts,rIphi_x,1.0);



	  // Loop over each quadrature point
	  for (int m=1; m<=mtmp; m++)
            {
	      // grid point x
	      xptsl.set( m, xc1 + 0.5*dxl*spts.get(m) );
              xptsr.set( m, xc2 + 0.5*dxr*spts.get(m) );

              
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

	      // Auxiliary values (aux) at each grid point
	      for (int ma=1; ma<=maux; ma++)
                {
		  auxvalsl.set(m,ma, 0.0 );

		  auxvalsr.set(m,ma, 0.0 );

		  for (int k=1; k<=mpoints; k++)
                    {
		      auxvalsl.set(m,ma, auxvalsl.get(m,ma) 
				  + phi.get(m,k) * auxI.get(1,iint,ma,k) );

		      auxvalsr.set(m,ma, auxvalsr.get(m,ma) 
				  + phi.get(m,k) * auxI.get(2,iint,ma,k) );
                    }
                }
            //printf("xcl=%e xtryl= %e \n",xc1,dxl);
            //printf("xcr=%e xtryr= %e \n",xc2,dxr);
            }

	  // Call user-supplied function to set fvals
          
	  Func(xptsl,qvalsl,auxvalsl,fvalsl);
          Func(xptsr,qvalsr,auxvalsr,fvalsr);
          /*
          for (int m=1; m<=mtmp; m++)
            {
               printf("xtryl= %e \n",xptsl.get(m));
               printf("xtryr= %e \n",xptsr.get(m));
               printf("qtryl= %e \n",qvalsl.get(m,1));
               printf("qtryr= %e \n",qvalsr.get(m,1));
            }*/


          //printf("i=%d iint=%d q= %e %e \n",i,iint,qI.get(1,iint,1,1),qI.get(1,iint,1,2)); 
          //printf("2i=%d iint=%d q= %e %e \n",i,iint,Iq.qget(2,iint,1,1),qI.get(2,iint,1,2)); 

          /*
	  for (inti m=1; m<=mtmp; m++)
            {
                printf("Points HERE! %d %e %e %e %e \n",m,xptsl.get(m),fvalsl.get(m,1),xptsr.get(m),fvalsl.get(m,1));
              
             }*/
 

	  // Evaluate integrals
	  if (mopt==0) // project onto Legendre basis
            {
	      for (int m1=1; m1<=mlength; m1++)
		for (int m2=1; m2<=mpoints; m2++)
		  {
		    double tmpl = 0.0;
                    double tmpr = 0.0;
		    for (int k=1; k<=mtmp; k++)
		      {
			tmpl += wgt.get(k)*fvalsl.get(k,m1)*phi.get(k,m2);
			tmpr += wgt.get(k)*fvalsr.get(k,m1)*phi.get(k,m2);
		      }
		    FI.set(1,iint,m1,m2, 0.5*tmpl );
                    FI.set(2,iint,m1,m2, 0.5*tmpr );
		  }
            }
	  else // project onto derivatives of Legendre basis
            {

              double Ul=auxI.get(1,iint,1,1);
              double Ur=auxI.get(2,iint,1,1);
              /*
	      for (int m1=1; m1<=mlength; m1++)             
		for (int m2=1; m2<=mpoints; m2++)
		  {
		    double tmpl = 0.0;
		    double tmpr = 0.0;
		    for (int k=1; k<=mtmp; k++)
		      {
			tmpl += wgt.get(k)*fvalsl.get(k,m1)*lIphi_x.get(k,m2);
			tmpr += wgt.get(k)*fvalsr.get(k,m1)*rIphi_x.get(k,m2);
		      }
		    FI.set(1,iint,m1,m2, 0.5*tmpl );
		    FI.set(2,iint,m1,m2, 0.5*tmpr );  
		  }*/

                for (int m2=1; m2<=mpoints; m2++)
                  for (int m3=1; m3<=mpoints; m3++)
                  {

                    double tmponl = 0.0;
                    double tmponr = 0.0;

                    double tmpIl = 0.0;
                    double tmpIr = 0.0;
                    for (int k=1; k<=mtmp; k++)
                      {
                        tmpIl += wgt.get(k)*Ul*phi.get(k,m3)*lIphi_x.get(k,m2);
                        tmpIr += wgt.get(k)*Ur*phi.get(k,m3)*rIphi_x.get(k,m2);
                        tmponl += wgt.get(k)*Ul*phi.get(k,m3)*phi_x.get(k,m2);
                        tmponr += wgt.get(k)*Ur*phi.get(k,m3)*phi_x.get(k,m2);
                      }
                    // Implicit.set(iint,m1,m2,m3, Implicit.get(iint,m1,m2,m3)-0.5*tmpIl );
                    //Implicit.set(iint,m1,kmax+m2,kmax+m3, Implicit.get(iint,m1,kmax+m2,kmax+m3)-0.5*tmpIr );
                    //Implicit.set(iint,m1,m2,m3, Implicit.get(iint,m1,m2,m3)+0.5*tmpIl );
                    //Implicit.set(iint,m1,kmax+m2,kmax+m3, Implicit.get(iint,m1,kmax+m2,kmax+m3)+0.5*tmpIr );
                    Implicit.set(iint,1,m2,m3, Implicit.get(iint,1,m2,m3)-0.5*tmpIl );
                    Implicit.set(iint,1,kmax+m2,kmax+m3, Implicit.get(iint,1,kmax+m2,kmax+m3)-0.5*tmpIl );
                    Implicit.set(iint,1,2*kmax+m2,2*kmax+m3, Implicit.get(iint,1,2*kmax+m2,2*kmax+m3)-0.5*tmpIl );
                    Implicit.set(iint,1,3*kmax+m2,3*kmax+m3, Implicit.get(iint,1,3*kmax+m2,3*kmax+m3)-0.5*tmpIr );
                    Implicit.set(iint,1,4*kmax+m2,4*kmax+m3, Implicit.get(iint,1,4*kmax+m2,4*kmax+m3)-0.5*tmpIr );
                    Implicit.set(iint,1,5*kmax+m2,5*kmax+m3, Implicit.get(iint,1,5*kmax+m2,5*kmax+m3)-0.5*tmpIr );
                    //if(abs(tmpIl)>1.0e-12 || abs(tmpIr)>1.0e-12)
                    //{printf("HERE2!!!! %e %e \n",0.5*tmpIl,0.5*tmpIr);}
                  }
                 /*
                 for(int m2=1;m2<=mpoints;m2++)
                 {
                  double tmp1=0.0;
                  double tmp2=0.0;
                  for (int m3=1;m3<=mpoints;m3++)
                  {
                    double tmpIl = 0.0;
                    double tmpIr = 0.0;
                    for (int k=1; k<=mtmp; k++)
                    {
                        tmpIl += wgt.get(k)*Ul*phi.get(k,m3)*lIphi_x.get(k,m2);
                        tmpIr += wgt.get(k)*Ur*phi.get(k,m3)*rIphi_x.get(k,m2);
                    }
                   tmp1=tmp1+0.5*tmpIl*qI.get(1,iint,1,m3);
                   tmp2=tmp2+0.5*tmpIr*qI.get(2,iint,1,m3);
                  }
                  //printf("HERE left! %e \n",tmp1-FI.get(1,1,1,m2));
                  //printf("HERE right! %e \n",tmp2-FI.get(2,1,1,m2));
                  }*/ 
            }
        }

      }
    }

}

void SampleBasisGrad_variable(const dTensor1& spts,
		     dTensor2& phi_x,double dx)
{
  const double odx = 1.0/dx;
  const int mpoints = phi_x.getsize(1);
  const int    kmax = phi_x.getsize(2);

  switch( kmax)
    {
    case 1:
      for (int m=1; m<=mpoints; m++)
	{
	  phi_x.set( m, 1, 0.0 );
	}
      break;

    case 2:
      for (int m=1; m<=mpoints; m++)
	{
	  phi_x.set( m, 1, 0.0 );
	  phi_x.set( m, 2, 2.0*sq3*odx );
	}
      break;

    case 3:
      for (int m=1; m<=mpoints; m++)
	{
	  const double xi  = spts.get(m);

	  phi_x.set( m, 1, 0.0 );
	  phi_x.set( m, 2, 2.0*sq3*odx );
	  phi_x.set( m, 3, 6.0*sq5*xi*odx );
	}
      break;

    case 4:
      for (int m=1; m<=mpoints; m++)
	{
	  const double xi  = spts.get(m);
	  const double xi2 = xi*xi;

	  phi_x.set( m, 1, 0.0 );
	  phi_x.set( m, 2, 2.0*sq3*odx );
	  phi_x.set( m, 3, 6.0*sq5*xi*odx );
	  phi_x.set( m, 4, 3.0*sq7*(5.0*xi2-1.0)*odx );
	}
      break;

    case 5:
      for (int m=1; m<=mpoints; m++)
	{
	  const double xi  = spts.get(m);
	  const double xi2 = xi*xi;
	  const double xi3 = xi2*xi;

	  phi_x.set( m, 1, 0.0 );
	  phi_x.set( m, 2, 2.0*sq3*odx );
	  phi_x.set( m, 3, 6.0*sq5*xi*odx );
	  phi_x.set( m, 4, 3.0*sq7*(5.0*xi2-1.0)*odx );
	  phi_x.set( m, 5, 15.0*xi* (7.0*xi2-3.0)*odx );
	}
      break;

    case 6:
      for (int m=1; m<=mpoints; m++)
	{
	  const double xi  = spts.get(m);
	  const double xi2 = xi*xi;
	  const double xi3 = xi2*xi;
	  const double xi4 = xi3*xi;

	  phi_x.set( m, 1, 0.0 );
	  phi_x.set( m, 2, 2.0*sq3*odx );
	  phi_x.set( m, 3, 6.0*sq5*xi*odx );
	  phi_x.set( m, 4, 3.0*sq7*(5.0*xi2-1.0)*odx );
	  phi_x.set( m, 5, 15.0*xi* (7.0*xi2-3.0)*odx );
	  phi_x.set( m, 6, (2.0*odx)*7.875*sq11*(5.0*xi4-10.0*onethird*xi2+5.0*onethird*oneseventh) );
	}
      break;

    }
}



