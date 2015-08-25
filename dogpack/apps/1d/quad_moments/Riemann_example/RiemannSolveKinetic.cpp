#include <cmath>
#include "dog_math.h"
#include "dogdefs.h"

// --------------------------------------------------------
// Kinetic Solver for the Bi-Gaussian Moment-Closure System
// --------------------------------------------------------
//
double RiemannSolve(const dTensor1& xedge,
                    const dTensor1& Ql,
                    const dTensor1& Qr,
                    const dTensor1& Auxl,
                    const dTensor1& Auxr,
                    dTensor1& Fl,
                    dTensor1& Fr,
                    void (*FluxFunc)(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&),
                    void (*SetWaveSpd)(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
                                       const dTensor1&,double&,double&))
{
  const int meqn = Ql.getsize();
  const int maux = Auxl.getsize();
  double smax_edge = 0.0e0;
  assert(meqn==5);
    
  // Calculate minimum and maximum speeds
  double s1,s2;
  SetWaveSpd(xedge,Ql,Qr,Auxl,Auxr,s1,s2);
  smax_edge = Max(fabs(s1),fabs(s2));

  // Calculate (+) and (-) fluxes using kinetic scheme
  dTensor1 Fp(meqn);
  void GetFluxPlus(const dTensor1& Ql, dTensor1& Fp);
  GetFluxPlus(Ql,Fp);

  dTensor1 Fm(meqn);
  void GetFluxMinus(const dTensor1& Qr, dTensor1& Fm);
  GetFluxMinus(Qr,Fm);
  

  // Calculate fluxes
  for (int m=1; m<=meqn; m++)
    {
      Fl.set(m, (Fp.get(m)+Fm.get(m)) );
      Fr.set(m, (Fp.get(m)+Fm.get(m)) );
    }
  
  return smax_edge;
}


void GetFluxPlus(const dTensor1& Q, dTensor1& Fp)
{
  dTensor1 BG(5);
  double MomentInversion(const dTensor1& Q, dTensor1& BG);
  const double alpha = MomentInversion(Q,BG);
  const double om1 = BG.get(1);
  const double om2 = BG.get(2);
  const double mu1 = BG.get(3);
  const double mu2 = BG.get(4);
  const double sig = BG.get(5);
  
  const double mu1_2 = mu1*mu1;
  const double mu1_3 = mu1_2*mu1;
  const double mu1_4 = mu1_3*mu1;
  const double mu1_5 = mu1_4*mu1;

  const double mu2_2 = mu2*mu2;
  const double mu2_3 = mu2_2*mu2;
  const double mu2_4 = mu2_3*mu2;
  const double mu2_5 = mu2_4*mu2;

  double S1p,S2p,S3p,S4p,S5p;
  if (sig<=1.0e-10)
    {
      double mu1p = Max(mu1,0.0);
      double mu2p = Max(mu2,0.0);

      S1p = om1*mu1p       + om2*mu2p;
      S2p = om1*mu1p*mu1   + om2*mu2p*mu2;
      S3p = om1*mu1p*mu1_2 + om2*mu2p*mu2_2;
      S4p = om1*mu1p*mu1_3 + om2*mu2p*mu2_3;
      S5p = om1*mu1p*mu1_4 + om2*mu2p*mu2_4;    
    }
  else
    {
      const double sig2  = sig*sig;
      
      const double Z1 = erf(mu1/sqrt(2.0*sig));
      const double Z2 = erf(mu2/sqrt(2.0*sig));
      const double G1 = exp(-mu1_2/(2.0*sig));
      const double G2 = exp(-mu2_2/(2.0*sig));
      
      S1p = (om1*(mu1 + sqrt(2.0*sig/pi)*G1 + mu1*Z1) +
	     om2*(mu2 + sqrt(2.0*sig/pi)*G2 + mu2*Z2))/2.0;
      
      S2p = (om1*(mu1_2 + sig + mu1*sqrt(2.0*sig/pi)*G1 + (mu1_2 + sig)*Z1) +
	     om2*(mu2_2 + sig + mu2*sqrt(2.0*sig/pi)*G2 + (mu2_2 + sig)*Z2))/2.0;
      
      S3p = (om1*(mu1_3 + 3.0*mu1*sig + (mu1_2 + 2.0*sig)*sqrt(2.0*sig/pi)*G1 + mu1*(mu1_2 + 3.0*sig)*Z1) +
	     om2*(mu2_3 + 3.0*mu2*sig + (mu2_2 + 2.0*sig)*sqrt(2.0*sig/pi)*G2 + mu2*(mu2_2 + 3.0*sig)*Z2))/2.0;
      
      S4p = (om1*(mu1_4 + 6.0*mu1_2*sig + 3.0*sig2 + mu1*(mu1_2 + 5.0*sig)*sqrt(2.0*sig/pi)*G1
		  + (mu1_4 + 6.0*sig*mu1_2 + 3.0*sig2)*Z1) +
	     om2*(mu2_4 + 6.0*mu2_2*sig + 3.0*sig2 + mu2*(mu2_2 + 5.0*sig)*sqrt(2.0*sig/pi)*G2
		  + (mu2_4 + 6.0*sig*mu2_2 + 3.0*sig2)*Z2))/2.0;
      
      S5p = (om1*(mu1_5 + 10.0*sig*mu1_3 + 15.0*mu1*sig2 +
		  (mu1_4 + 9.0*mu1_2*sig + 8.0*sig2)*sqrt(2.0*sig/pi)*G1 +
		  mu1*(mu1_4 + 10.0*mu1_2*sig + 15.0*sig2)*Z1) +
	     om2*(mu2_5 + 10.0*sig*mu2_3 + 15.0*mu2*sig2 +
		  (mu2_4 + 9.0*mu2_2*sig + 8.0*sig2)*sqrt(2.0*sig/pi)*G2 +
		  mu2*(mu2_4 + 10.0*mu2_2*sig + 15.0*sig2)*Z2))/2.0;
    }
      
  Fp.set(1, S1p );
  Fp.set(2, S2p );
  Fp.set(3, S3p );
  Fp.set(4, S4p );
  Fp.set(5, S5p );
}


void GetFluxMinus(const dTensor1& Q, dTensor1& Fm)
{
  dTensor1 BG(5);
  double MomentInversion(const dTensor1& Q, dTensor1& BG);
  const double alpha = MomentInversion(Q,BG);
  const double om1 = BG.get(1);
  const double om2 = BG.get(2);
  const double mu1 = BG.get(3);
  const double mu2 = BG.get(4);
  const double sig = BG.get(5);

  const double mu1_2 = mu1*mu1;
  const double mu1_3 = mu1_2*mu1;
  const double mu1_4 = mu1_3*mu1;
  const double mu1_5 = mu1_4*mu1;

  const double mu2_2 = mu2*mu2;
  const double mu2_3 = mu2_2*mu2;
  const double mu2_4 = mu2_3*mu2;
  const double mu2_5 = mu2_4*mu2;

  double S1m,S2m,S3m,S4m,S5m;
  if (sig<=1.0e-10)
    {
      double mu1m = Min(mu1,0.0);
      double mu2m = Min(mu2,0.0);

      S1m = om1*mu1m       + om2*mu2m;
      S2m = om1*mu1m*mu1   + om2*mu2m*mu2;
      S3m = om1*mu1m*mu1_2 + om2*mu2m*mu2_2;
      S4m = om1*mu1m*mu1_3 + om2*mu2m*mu2_3;
      S5m = om1*mu1m*mu1_4 + om2*mu2m*mu2_4;    
    }
  else
    {
      const double sig2  = sig*sig;
      
      const double Z1 = -erf(mu1/sqrt(2.0*sig));
      const double Z2 = -erf(mu2/sqrt(2.0*sig));
      const double G1 = -exp(-mu1_2/(2.0*sig));
      const double G2 = -exp(-mu2_2/(2.0*sig));
      
      S1m = (om1*(mu1 + sqrt(2.0*sig/pi)*G1 + mu1*Z1) +
	     om2*(mu2 + sqrt(2.0*sig/pi)*G2 + mu2*Z2))/2.0;
      
      S2m = (om1*(mu1_2 + sig + mu1*sqrt(2.0*sig/pi)*G1 + (mu1_2 + sig)*Z1) +
	     om2*(mu2_2 + sig + mu2*sqrt(2.0*sig/pi)*G2 + (mu2_2 + sig)*Z2))/2.0;
      
      S3m = (om1*(mu1_3 + 3.0*mu1*sig + (mu1_2 + 2.0*sig)*sqrt(2.0*sig/pi)*G1 + mu1*(mu1_2 + 3.0*sig)*Z1) +
	     om2*(mu2_3 + 3.0*mu2*sig + (mu2_2 + 2.0*sig)*sqrt(2.0*sig/pi)*G2 + mu2*(mu2_2 + 3.0*sig)*Z2))/2.0;
      
      S4m = (om1*(mu1_4 + 6.0*mu1_2*sig + 3.0*sig2 + mu1*(mu1_2 + 5.0*sig)*sqrt(2.0*sig/pi)*G1
		  + (mu1_4 + 6.0*sig*mu1_2 + 3.0*sig2)*Z1) +
	     om2*(mu2_4 + 6.0*mu2_2*sig + 3.0*sig2 + mu2*(mu2_2 + 5.0*sig)*sqrt(2.0*sig/pi)*G2
		  + (mu2_4 + 6.0*sig*mu2_2 + 3.0*sig2)*Z2))/2.0;
      
      S5m = (om1*(mu1_5 + 10.0*sig*mu1_3 + 15.0*mu1*sig2 +
		  (mu1_4 + 9.0*mu1_2*sig + 8.0*sig2)*sqrt(2.0*sig/pi)*G1 +
		  mu1*(mu1_4 + 10.0*mu1_2*sig + 15.0*sig2)*Z1) +
	     om2*(mu2_5 + 10.0*sig*mu2_3 + 15.0*mu2*sig2 +
		  (mu2_4 + 9.0*mu2_2*sig + 8.0*sig2)*sqrt(2.0*sig/pi)*G2 +
		  mu2*(mu2_4 + 10.0*mu2_2*sig + 15.0*sig2)*Z2))/2.0;
    }

  Fm.set(1, S1m );
  Fm.set(2, S2m );
  Fm.set(3, S3m );
  Fm.set(4, S4m );
  Fm.set(5, S5m );

}
