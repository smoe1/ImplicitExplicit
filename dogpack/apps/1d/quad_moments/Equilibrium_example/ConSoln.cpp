#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include "dogdefs.h"
#include "DogState1d.h"
#include "QuadMomentParams.h" 
using namespace std;

void ConSoln(const int method[], 
	     const dTensor2& node, 
	     const dTensorBC3& aux,
	     const dTensorBC3& q, 
	     double t, 
	     string outputdir)
{
  const int melems = q.getsize(1);
  const int meqn   = q.getsize(2);
  const int kmax   = q.getsize(3);
  const int mbc    = q.getmbc();
  const int maux   = aux.getsize(2);
  string fname1 = outputdir+"/conservation.dat";
  string fname2 = outputdir+"/equilibrium.dat";
  ofstream write_file1,write_file2;
  dTensor1 qsum(meqn);
  dTensor1 res_sum(meqn);
  
  if (t==0) 
    {
      write_file1.open(fname1.c_str(), ofstream::out);
      write_file2.open(fname2.c_str(), ofstream::out);

      const double epsilon = quadMomentParams.epsilon;
      write_file2 << setprecision(16);
      write_file2 << setw(24) << scientific << epsilon << endl;
    }
  else
    {
      write_file1.open(fname1.c_str(), ofstream::app);
      write_file2.open(fname2.c_str(), ofstream::app);
    }
  
  // -----------------
  // CONSERVATION
  // -----------------
  for (int m=1; m<=meqn; m++)
    {
      qsum.set(m,0.0);
      
      for (int i=1; i<=melems; i++)
	{
	  double x = node.get(i,1);
	  double dtmp = node.get(i+1,1)-node.get(i,1);
	  double qtmp = q.get(i,m,1);
	  
	  qsum.set(m, (qsum.get(m) + dtmp*qtmp) );
	}
    }

  // ---------------------
  // EQUILIBRIUM CHECK
  // ---------------------
  dTensorBC3 Evals(melems, 1, kmax, mbc);
  dTensorBC3 auxtmp(melems, maux, kmax, mbc);
  for (int i=(1-mbc); i<=(melems+mbc); i++)
    for (int m=1; m<=maux; m++)
      for (int k=1; k<=kmax; k++)
	{
	  auxtmp.set(i,m,k, aux.get(i,m,k) );
	}
  void ComputeElecField(double t, 
			const dTensor2& node,
			const dTensorBC3& qvals, 
			dTensorBC3& auxtmp, 
			dTensorBC3& Evals);
  ComputeElecField(fetch_dogState().get_time(), node, q, auxtmp, Evals);

  double M1diff = 0.0;
  for (int i=1; i<=melems; i++)
    {
      double rho1 = q.get(i,1,1);
      double rho2 = q.get(i,1,2);
      double rho3 = q.get(i,1,3);

      double M1   = q.get(i,2,1);
      double M2   = q.get(i,2,2);
      double M3   = q.get(i,2,3);

      double E1   = -Evals.get(i,1,1);
      double E2   = -Evals.get(i,1,2);
      double E3   = -Evals.get(i,1,3);

      double rhoE1 = rho1*E1 + rho2*E2 + rho3*E3;
      double rhoE2 = rho2*E1 + ((2.0/5.0)*rho3*sq5+rho1)*E2 
	+ (2.0/5.0)*E3*sq5*rho2;
      double rhoE3 = (2.0/7.0)*E3*sq5*rho3 + rho3*E1 
	+ (2.0/5.0)*sq5*E2*rho2 + E3*rho1;
      
      M1diff = M1diff 
	+ pow(M1+rhoE1,2) + pow(M2+rhoE2,2) + pow(M3+rhoE3,2);
    }
  double dx = node.get(2,1)-node.get(1,1);
  M1diff = sqrt(dx*M1diff);

  write_file1 << setprecision(16);
  write_file1 << setw(24) << scientific << t << " ";
  for (int m=1; m<=meqn; m++)
    {
      if (abs(qsum.get(m)) < 1.0e-99) {qsum.set(m, 0.0);}
      write_file1 << setw(24) << scientific << qsum.get(m) << " ";
    }
  write_file1 << endl;

  write_file2 << setprecision(16);
  write_file2 << setw(24) << scientific << t << " ";
  if (fabs(M1diff) < 1.0e-99) 
    {M1diff = 0.0;}
  write_file2 << setw(24) << scientific << M1diff << " ";
  write_file2 << endl;
  
}
