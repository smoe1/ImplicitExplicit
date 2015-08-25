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
  assert(meqn==4);
    
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
  void MomentInversion(const dTensor1& Q, dTensor1& BG);
  MomentInversion(Q,BG);
  const double om1 = BG.get(1);
  const double om2 = BG.get(2);
  const double mu1 = BG.get(3);
  const double mu2 = BG.get(4);
  const double sig = BG.get(5);
  const double ssig = sqrt(sig);
  const double ss2 = 0.5*ssig;

  const double mu1_2 = mu1*mu1;
  const double mu2_2 = mu2*mu2;
  const double mu1_3 = mu1*mu1_2;
  const double mu2_3 = mu2*mu2_2;
  const double mu1_4 = mu1*mu1_3;
  const double mu2_4 = mu2*mu2_3;

  double S1p,S2p,S3p,S4p;
  double T1p,T2p,T3p,T4p;

  if (mu1-ss2>=0.0)
    {
      S1p = om1;
      S2p = om1*mu1;
      S3p = om1*mu1_2+om1*sig/24.0;
      S4p = om1*mu1_3+om1*mu1*sig/8.0;
    }
  else if (mu1>=0.0)
    {
      S1p = 2.0*mu1*om1*(ssig-mu1)/sig 
	+ 0.5*om1;
      S2p = onethird*om1*mu1_2*(3.0*ssig-2.0*mu1)/sig
	+ 0.25*onethird*om1*(6.0*mu1+ssig);
      S3p = onethird*om1*mu1_3*(2.0*ssig-mu1)/sig 
	+ 0.0625*onethird*om1*(24.0*mu1_2+8.0*mu1*ssig+sig);
      S4p = 0.1*om1*mu1_4*(5.0*ssig-2.0*mu1)/sig 
	+ 0.00625*om1*(80.0*mu1_3+40.0*mu1_2*ssig+10.0*sig*mu1+pow(sig,1.5));
    }
  else if (mu1+ss2>=0.0)
    {
      S1p = 0.5*pow(2.0*mu1+ssig,2)*om1/sig;
      S2p = 0.25*onethird*pow(2.0*mu1+ssig,3)*om1/sig;
      S3p = 0.0625*onethird*om1*pow(2.0*mu1+ssig,4)/sig;
      S4p = 0.00625*om1*pow(2.0*mu1+ssig,5)/sig;
    }
  else
    {
      S1p = 0.0;
      S2p = 0.0;
      S3p = 0.0;
      S4p = 0.0;
    }

  if (mu2-ss2>=0.0)
    {
      T1p = om2;
      T2p = om2*mu2;
      T3p = om2*mu2_2+om2*sig/24.0;
      T4p = om2*mu2_3+om2*mu2*sig/8.0;
    }
  else if (mu2>=0.0)
    {
      T1p = 2.0*mu2*om2*(ssig-mu2)/sig 
	+ 0.5*om2;
      T2p = onethird*om2*mu2_2*(3.0*ssig-2.0*mu2)/sig
	+ 0.25*onethird*om2*(6.0*mu2+ssig);
      T3p = onethird*om2*mu2_3*(2.0*ssig-mu2)/sig 
	+ 0.0625*onethird*om2*(24.0*mu2_2+8.0*mu2*ssig+sig);
      T4p = 0.1*om2*mu2_4*(5.0*ssig-2.0*mu2)/sig 
	+ 0.00625*om2*(80.0*mu2_3+40.0*mu2_2*ssig+10.0*sig*mu2+pow(sig,1.5));
    }
  else if (mu2+ss2>=0.0)
    {
      T1p = 0.5*pow(2.0*mu2+ssig,2)*om2/sig;
      T2p = 0.25*onethird*pow(2.0*mu2+ssig,3)*om2/sig;
      T3p = 0.0625*onethird*om2*pow(2.0*mu2+ssig,4)/sig;
      T4p = 0.00625*om2*pow(2.0*mu2+ssig,5)/sig;
    }
  else
    {
      T1p = 0.0;
      T2p = 0.0;
      T3p = 0.0;
      T4p = 0.0;
    }

  /*
  printf("  om1 = %e\n",om1);
  printf("  mu1 = %e\n",mu1);
  printf("  om2 = %e\n",om2);
  printf("  mu2 = %e\n",mu2);
  printf("  sig = %e\n",sig);
  printf("\n");
  printf("  S1p = %e\n",S1p);
  printf("  S2p = %e\n",S2p);
  printf("  S3p = %e\n",S3p);
  printf("  S4p = %e\n",S4p);
  printf("\n");
  printf("  T1p = %e\n",T1p);
  printf("  T2p = %e\n",T2p);
  printf("  T3p = %e\n",T3p);
  printf("  T4p = %e\n",T4p);
  printf("\n");
  */

  Fp.set(1, S1p+T1p );
  Fp.set(2, S2p+T2p );
  Fp.set(3, S3p+T3p );
  Fp.set(4, S4p+T4p );
}


void GetFluxMinus(const dTensor1& Q, dTensor1& Fm)
{
  dTensor1 BG(5);
  void MomentInversion(const dTensor1& Q, dTensor1& BG);
  MomentInversion(Q,BG);
  const double om1 = BG.get(1);
  const double om2 = BG.get(2);
  const double mu1 = BG.get(3);
  const double mu2 = BG.get(4);
  const double sig = BG.get(5);
  const double ssig = sqrt(sig);
  const double ss2 = 0.5*ssig;

  const double mu1_2 = mu1*mu1;
  const double mu2_2 = mu2*mu2;
  const double mu1_3 = mu1*mu1_2;
  const double mu2_3 = mu2*mu2_2;
  const double mu1_4 = mu1*mu1_3;
  const double mu2_4 = mu2*mu2_3;

  double S1m,S2m,S3m,S4m;
  double T1m,T2m,T3m,T4m;

  if (mu1+ss2<=0.0)
    {      
      S1m = om1;
      S2m = om1*mu1;
      S3m = om1*mu1_2+om1*sig/24.0;
      S4m = om1*mu1_3+om1*mu1*sig/8.0;
    }
  else if (mu1<=0.0)
    {      
      S1m = 0.5*om1
	- 2.0*mu1*om1*(mu1+ssig)/sig;
      S2m = -0.25*onethird*om1*(-6.0*mu1+ssig)
	- onethird*om1*mu1_2*(2.0*mu1+3.0*ssig)/sig;
      S3m = -0.0625*onethird*om1*(-24.0*mu1_2+8.0*mu1*ssig-sig)
	- onethird*om1*mu1_3*(mu1+2.0*ssig)/sig;
      S4m = -0.00625*om1*(-80.0*mu1_3+40.0*mu1_2*ssig-10.0*sig*mu1+pow(sig,1.5))
	-0.1*om1*mu1_4*(2.0*mu1+5*ssig)/sig;
    }
  else if (mu1-ss2<=0.0)
    {      
      S1m = 0.5*pow(ssig-2.0*mu1,2)*om1/sig;
      S2m = -0.25*onethird*pow(ssig-2.0*mu1,3)*om1/sig;
      S3m = 0.0625*onethird*om1*pow(ssig-2.0*mu1,4)/sig;
      S4m = -0.00625*om1*pow(ssig-2.0*mu1,5)/sig;
    }
  else
    {
      S1m = 0.0;
      S2m = 0.0;
      S3m = 0.0;
      S4m = 0.0;
    }

  if (mu2+ss2<=0.0)
    {
      T1m = om2;
      T2m = om2*mu2;
      T3m = om2*mu2_2+om2*sig/24.0;
      T4m = om2*mu2_3+om2*mu2*sig/8.0;
    }
  else if (mu2<=0.0)
    {
      T1m = 0.5*om2
	- 2.0*mu2*om2*(mu2+ssig)/sig;
      T2m = -0.25*onethird*om2*(-6.0*mu2+ssig)
	- onethird*om2*mu2_2*(2.0*mu2+3.0*ssig)/sig;
      T3m = -0.0625*onethird*om2*(-24.0*mu2_2+8.0*mu2*ssig-sig)
	- onethird*om2*mu2_3*(mu2+2.0*ssig)/sig;
      T4m = -0.00625*om2*(-80.0*mu2_3+40.0*mu2_2*ssig-10.0*sig*mu2+pow(sig,1.5))
	-0.1*om2*mu2_4*(2.0*mu2+5*ssig)/sig;
    }
  else if (mu2-ss2<=0.0)
    {
      T1m = 0.5*pow(ssig-2.0*mu2,2)*om2/sig;
      T2m = -0.25*onethird*pow(ssig-2.0*mu2,3)*om2/sig;
      T3m = 0.0625*onethird*om2*pow(ssig-2.0*mu2,4)/sig;
      T4m = -0.00625*om2*pow(ssig-2.0*mu2,5)/sig;
    }
  else
    {
      T1m = 0.0;
      T2m = 0.0;
      T3m = 0.0;
      T4m = 0.0;
    }  

  /*
  printf("  om1 = %e\n",om1);
  printf("  mu1 = %e\n",mu1);
  printf("  om2 = %e\n",om2);
  printf("  mu2 = %e\n",mu2);
  printf("  sig = %e\n",sig);
  printf("\n");
  printf("  S1m = %e\n",S1m);
  printf("  S2m = %e\n",S2m);
  printf("  S3m = %e\n",S3m);
  printf("  S4m = %e\n",S4m);
  printf("\n");
  printf("  T1m = %e\n",T1m);
  printf("  T2m = %e\n",T2m);
  printf("  T3m = %e\n",T3m);
  printf("  T4m = %e\n",T4m);
  printf("\n");
  */

  Fm.set(1, S1m+T1m );
  Fm.set(2, S2m+T2m );
  Fm.set(3, S3m+T3m );
  Fm.set(4, S4m+T4m );
}
