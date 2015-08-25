#include <cmath>
#include "dog_math.h"
#include "tensors.h"
#include "dogdefs.h"


// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
//
void SetWaveSpd(const dTensor1& xedge,
		const dTensor1& Ql,
		const dTensor1& Qr,
		const dTensor1& Auxl,
               const dTensor1& Auxr,
		double& s1,double& s2)
{    
      void LinTrans(const double* Q, const double* L, double* Qt);
      void GetMus(const double rho, const double* u,
                  const double* L, double* Qt,
                  double mut[][3], double mu[][3]);
   
  // left wave speed
  /*
      const double M000l = Ql.get(1);
      const double M100l = Ql.get(2);
      const double M010l = Ql.get(3);
      const double M001l = Ql.get(4);
      const double M200l = Ql.get(5);
      const double M110l = Ql.get(6);
      const double M101l = Ql.get(7);
      const double M020l = Ql.get(8);
      const double M011l = Ql.get(9);
      const double M002l = Ql.get(10);
      const double M300l = Ql.get(11);
      const double M210l = Ql.get(12);
      const double M201l = Ql.get(13);
      const double M120l = Ql.get(14);
      const double M111l = Ql.get(15);
      const double M102l = Ql.get(16);
      const double M030l = Ql.get(17);
      const double M021l = Ql.get(18);
      const double M012l = Ql.get(19);
      const double M003l = Ql.get(20);  */

      const double M000l = Qr.get(1);
      const double M100l = Qr.get(2);
      const double M010l = Qr.get(3);
      const double M001l = Qr.get(4);
      const double M200l = Qr.get(5);
      const double M110l = Qr.get(6);
      const double M101l = Qr.get(7);
      const double M020l = Qr.get(8);
      const double M011l = Qr.get(9);
      const double M002l = Qr.get(10);
      const double M300l = Qr.get(11);
      const double M030l = Qr.get(12);
      const double M003l = Qr.get(13);
      const double M111l = Qr.get(14);
  

      const double rhol  = M000l;
      const double u1l   = M100l/rhol;
      const double u2l   = M010l/rhol;
      const double u3l   = M001l/rhol;
      const double P11l  = M200l - rhol*u1l*u1l;
      const double P22l  = M020l - rhol*u2l*u2l;
      const double P33l  = M002l - rhol*u3l*u3l;
      const double P12l  = M110l - rhol*u1l*u2l;
      const double P13l  = M101l - rhol*u1l*u3l;
      const double P23l  = M011l - rhol*u2l*u3l;
      const double Q111l = M300l - 3.0*u1l*P11l - rhol*pow(u1l,3);
      //const double Q112l = M210l - 2.0*u1l*P12l - u2l*P11l - rhol*u1l*u1l*u2l;
      //const double Q113l = M201l - 2.0*u1l*P13l - u3l*P11l - rhol*u1l*u1l*u3l;
      //const double Q122l = M120l - 2.0*u2l*P12l - u1l*P22l - rhol*u1l*u2l*u2l;
      const double Q123l = M111l - u1l*P23l - u2l*P13l - u3l*P12l - rhol*u1l*u2l*u3l;
      //const double Q133l = M102l - 2.0*u3l*P13l - u1l*P33l - rhol*u1l*u3l*u3l;
      const double Q222l = M030l - 3.0*u2l*P22l - rhol*pow(u2l,3);
      //const double Q223l = M021l - 2.0*u2l*P23l - u3l*P22l - rhol*u2l*u2l*u3l;
      //const double Q233l = M012l - 2.0*u3l*P23l - u2l*P33l - rhol*u2l*u3l*u3l;
      const double Q333l = M003l - 3.0*u3l*P33l - rhol*pow(u3l,3);

      const double detPl= P11l*P22l*P33l-P11l*P23l*P23l-P22l*P13l*P13l-P33l*P12l*P12l+2.0*P12l*P13l*P23l;

      if ( rhol<=0.0 || P11l<0.0 || P22l<0.0 || P33l<0.0 || P11l*P22l-P12l*P12l<0.0 || detPl<0.0)
        {
          printf("              rhol = %22.16e \n",rhol);
          printf("              P11l = %22.16e \n",P11l);
          printf("              P22l = %22.16e \n",P22l);
          printf("              P33l = %22.16e \n",P33l);
          printf("  P11l*P22l-P12l*P12l = %22.16e \n",P11l*P22l-P12l*P12l);
          printf("             detPl = %22.16e \n",detPl);
          printf("\n");
          exit(1);
        }

      const double L11l = sqrt(P11l);
      const double L22l = sqrt(P22l-P12l*P12l/P11l);
      const double L33l = sqrt(P33l-P13l*P13l/P11l-pow((P11l*P23l-P13l*P13l),2)/(P11l*P11l*P22l-P11l*P12l*P12l));
      const double L21l = P12l/sqrt(P11l);
      const double L31l = P13l/sqrt(P11l);
      const double L32l = (P11l*P23l-P12l*P13l)/sqrt(P11l*P11l*P22l-P11l*P12l*P12l);

      const double ul[3]   = {u1l,u2l,u3l};
      const double Ll[6]   = {L11l,L22l,L33l,L21l,L31l,L32l};
      //const double Qll[10]  = {Q111l,Q112l,Q113l,Q122l,Q123l,Q133l,Q222l,Q223l,Q233l,Q333l};
      const double Qll[4]  = {Q111l,Q222l,Q333l,Q123l};
      //double Qtl[10];
      double Qtl[4];

      LinTrans(Qll,Ll,Qtl);
           
      double mutl[8][3];
      double mul[8][3];
      
      GetMus(rhol,ul,Ll,Qtl,mutl,mul);


      s1 = abs(mul[0][0]);

      for (int i=1;i<8;i++)
      {
        s1 = -Max(s1,abs(mul[i][0]));
      }

   // right wave speed
  /*
      const double M000r = Qr.get(1);
      const double M100r = Qr.get(2);
      const double M010r = Qr.get(3);
      const double M001r = Qr.get(4);
      const double M200r = Qr.get(5);
      const double M110r = Qr.get(6);
      const double M101r = Qr.get(7);
      const double M020r = Qr.get(8);
      const double M011r = Qr.get(9);
      const double M002r = Qr.get(10);
      const double M300r = Qr.get(11);
      const double M210r = Qr.get(12);
      const double M201r = Qr.get(13);
      const double M120r = Qr.get(14);
      const double M111r = Qr.get(15);
      const double M102r = Qr.get(16);
      const double M030r = Qr.get(17);
      const double M021r = Qr.get(18);
      const double M012r = Qr.get(19);
      const double M003r = Qr.get(20);  */

      const double M000r = Qr.get(1);
      const double M100r = Qr.get(2);
      const double M010r = Qr.get(3);
      const double M001r = Qr.get(4);
      const double M200r = Qr.get(5);
      const double M110r = Qr.get(6);
      const double M101r = Qr.get(7);
      const double M020r = Qr.get(8);
      const double M011r = Qr.get(9);
      const double M002r = Qr.get(10);
      const double M300r = Qr.get(11);
      const double M030r = Qr.get(12);
      const double M003r = Qr.get(13);
      const double M111r = Qr.get(14);

      const double rhor  = M000r;
      const double u1r   = M100r/rhor;
      const double u2r   = M010r/rhor;
      const double u3r   = M001r/rhor;
      const double P11r  = M200r - rhor*u1r*u1r;
      const double P22r  = M020r - rhor*u2r*u2r;
      const double P33r  = M002r - rhor*u3r*u3r;
      const double P12r  = M110r - rhor*u1r*u2r;
      const double P13r  = M101r - rhor*u1r*u3r;
      const double P23r  = M011r - rhor*u2r*u3r;
      const double Q111r = M300r - 3.0*u1r*P11r - rhor*pow(u1r,3);
      //const double Q112r = M210r - 2.0*u1r*P12r - u2r*P11r - rhor*u1r*u1r*u2r;
      //const double Q113r = M201r - 2.0*u1r*P13r - u3r*P11r - rhor*u1r*u1r*u3r;
      //const double Q122r = M120r - 2.0*u2r*P12r - u1r*P22r - rhor*u1r*u2r*u2r;
      const double Q123r = M111r - u1r*P23r - u2r*P13r - u3r*P12r - rhor*u1r*u2r*u3r;
      //const double Q133r = M102r - 2.0*u3r*P13r - u1r*P33r - rhor*u1r*u3r*u3r;
      const double Q222r = M030r - 3.0*u2r*P22r - rhor*pow(u2r,3);
      //const double Q223r = M021r - 2.0*u2r*P23r - u3l*P22r - rhor*u2r*u2r*u3r;
      //const double Q233r = M012r - 2.0*u3r*P23r - u2l*P33r - rhor*u2r*u3r*u3r;
      const double Q333r = M003r - 3.0*u3r*P33r - rhor*pow(u3r,3);

      const double detPr= P11r*P22r*P33r-P11r*P23r*P23r-P22r*P13r*P13r-P33r*P12r*P12r+2.0*P12r*P13r*P23r;

      if ( rhor<=0.0 || P11r<0.0 || P22r<0.0 || P33r<0.0 || P11r*P22r-P12r*P12r<0.0 || detPr<0.0)
        {
          printf("              rhor = %22.16e \n",rhor);
          printf("              P11r = %22.16e \n",P11r);
          printf("              P22r = %22.16e \n",P22r);
          printf("              P33r = %22.16e \n",P33r);
          printf("  P11r*P22r-P12r*P12r = %22.16e \n",P11r*P22r-P12r*P12r);
          printf("             detPr = %22.16e \n",detPr);
          printf("\n");
          exit(1);
        }

      const double L11r = sqrt(P11r);
      const double L22r = sqrt(P22r-P12r*P12r/P11r);
      const double L33r = sqrt(P33r-P13r*P13r/P11r-pow((P11r*P23r-P13r*P13r),2)/(P11r*P11r*P22r-P11r*P12r*P12r));
      const double L21r = P12r/sqrt(P11r);
      const double L31r = P13r/sqrt(P11r);
      const double L32r = (P11r*P23r-P12r*P13r)/sqrt(P11r*P11r*P22r-P11r*P12r*P12r);

      const double ur[3]   = {u1r,u2r,u3r};
      const double Lr[6]   = {L11r,L22r,L33r,L21r,L31r,L32r};
      //const double Qrr[10]   = {Q111r,Q112r,Q113r,Q122r,Q123r,Q133r,Q222r,Q223r,Q233r,Q333r};
      const double Qrr[4]  = {Q111r,Q222r,Q333r,Q123r};
      //double Qtr[10];
      double Qtr[4];

      LinTrans(Qrr,Lr,Qtr);
           
      double mutr[8][3];
      double mur[8][3];
      
      GetMus(rhor,ur,Lr,Qtr,mutr,mur); 


      s2 = abs(mur[0][0]);

      for (int i=1;i<8;i++)
      {
        s2 = Max(s2,abs(mur[i][0]));
      }
 
}



