#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "dogdefs.h"
using namespace std;

// =========================================================================
//
//  Copyright J.A. Rossmanith
//
//  This software is made available for research and instructional use only.
//  You may copy and use this software without charge for these non-commercial
//  purposes, provided that the copyright notice and associated text is
//  reproduced on all copies.  For all other uses (including distribution of
//  modified versions), please contact the author at the address given below.
//
//  *** This software is made available "as is" without any assurance that it
//  *** will work for your purposes.  The software may in fact have defects, so
//  *** use the software at your own risk.
//
//  -------------------------------
//  DoGPack
//  -------------------------------
//
//    Lead Developer:  
//             James Rossmanith
//             Iowa State University
//             Department of Mathematics
//             396 Carver Hall
//             Ames, IA 50011
//             rossmani@iastate.edu
// =========================================================================

int main()
{  
    // ------------------------------------------------------------
    // Function definitions
    double bottom(double);
    double Newton(double,double,double,double,double,int);
    // ------------------------------------------------------------

    int i;
    
    int numpts;
    int mcrit;
    int mstrength;
    double ShockLoc,ShockStr;
    double Froude,Qflow,Const;
    double a0,a1,a2,root;
    
    // --------------------------------------------
    // "param.data" is the input file for this code
    // --------------------------------------------
    ifstream read_file1("param.data", ios::in);
    read_file1 >> numpts;
    int nmid = (numpts+2)/2;
    double hmid = pow(0.25,2.0/3.0); //0.35;
    double x=0.0;
    double dx = 1.0/double(numpts);

    dTensor1 height(numpts+1);
    dTensor1 xout(numpts+1);
    dTensor1 b(numpts+1);
  
    read_file1 >> mcrit;  // if mcrit==0, then flow is symmetric,
                          // if mcrit!=0, then flow has a shock
    read_file1 >> ShockLoc;
    read_file1 >> mstrength;
    read_file1.close();

    // initialize x and b
    xout.set(1,x);
    b.set(1,bottom(x));
    for (i=1; i<=numpts; i++)
    {
        x = x + dx;
	xout.set(i+1,x);
	b.set(i+1,bottom(x));
    }

    // Mid-point
    height.set(nmid,hmid);

    // Compute important constants
    if (mcrit == 0)
    { 
        Froude = 0.8;
	Qflow = pow(hmid,1.5e0)*Froude;
	Const = hmid + b.get(nmid) + 0.5*(Qflow*Qflow)/(hmid*hmid); 
    }  
    else
    {
        Froude = 1.0;
	Qflow = pow(hmid,1.5e0)*Froude;
	Const = hmid + b.get(nmid) + 0.5*(Qflow*Qflow)/(hmid*hmid);
    }

    // Subcritical flow
    if (mcrit == 0)
    {
        // Loop in backward direction
        for (i=nmid-1; i>=1; i--)
	{
	    a0 = 0.5*Qflow*Qflow;
	    a1 = 0.0e0;
	    a2 = b.get(i)-Const;
	  
	    root = Newton(a0,a1,a2,height.get(i+1),1.0e-12,500);
	  
	    height.set(i,root);
	}
      
	// Loop in forward direction
	for (i=nmid+1; i<=(numpts+1); i++)
	{
	    a0 = 0.5*Qflow*Qflow;
	    a1 = 0.0e0;
	    a2 = b.get(i)-Const;
	  
	    root = Newton(a0,a1,a2,height.get(i-1),1.0e-12,500);
	  
	    height.set(i,root);
	}
    }

    // Transcritical flow
    if (mcrit == 1)
    {
        // Loop in backward direction (the subcritical direction)
        a0 = 0.5*Qflow*Qflow;
	a1 = 0.0e0;
	a2 = b.get(nmid-1)-Const;      
	root = Newton(a0,a1,a2,(1.1)*height.get(nmid),1.0e-12,500);	  
	height.set(nmid-1,root);

	for (i=nmid-2; i>=1; i--)
	{
	    a0 = 0.5*Qflow*Qflow;
	    a1 = 0.0e0;
	    a2 = b.get(i)-Const;
	  
	    root = Newton(a0,a1,a2,height.get(i+1),1.0e-12,500);
	  
	    height.set(i,root);
	}
      
	// Loop in forward direction up to shock (the supercritical regime)

	a0 = 0.5*Qflow*Qflow;
	a1 = 0.0e0;
	a2 = b.get(nmid+1)-Const;      
	root = Newton(a0,a1,a2,(0.9)*height.get(nmid),1.0e-12,500);	  
	height.set(nmid+1,root);

	i = nmid+1;
	x = xout.get(i);

	while( i<=(numpts+1) && x<ShockLoc )
	{
	    i = i+1;
	    x = xout.get(i);
	    a0 = 0.5*Qflow*Qflow;
	    a1 = 0.0e0;
	    a2 = b.get(i)-Const;
	    
	    root = Newton(a0,a1,a2,height.get(i-1),1.0e-12,500);
	  
	    height.set(i,root);
	}

	int mshock = i;
      
	// compute shock strength based on Rankine-Hugoniot conditions
	double hnew;
	if (mstrength==0)
	{ hnew = height.get(mshock); }
	else
	{
	    double Fr_tmp = 0.01;
	    double hguess = pow((Qflow/Fr_tmp), (2.0/3.0));
	    double tmp = (Qflow*Qflow)/(height.get(mshock)) 
	      + 0.5*height.get(mshock)*height.get(mshock);
	  
	    a0 = 2.0*Qflow*Qflow;
	    a1 = -2.0*tmp;
	    a2 = 0.0e0;
	    
	    hnew = Newton(a0,a1,a2,height.get(1),1.0e-12,500);

	}
	 
	Const = hnew + b.get(mshock) + 0.5*(Qflow*Qflow)/(hnew*hnew);
      
	height.set(mshock+1,hnew);
	for (i=mshock+2; i<=(numpts+1); i++)
	{
	    a0 = 0.5*Qflow*Qflow;
	    a1 = 0.0e0;
	    a2 = b.get(i)-Const;
	  
	    root = Newton(a0,a1,a2,height.get(i-1),1.0e-12,500);
	  
	    height.set(i,root);
	}
      
    }
    
    // Output exact solution
    ofstream write_file("shllw.dat", ios::out);   
    write_file << setprecision(16);
    for (i=1; i<=(numpts+1); i++)
    {
        write_file << setw(24) << scientific << xout.get(i) << " ";
	write_file << setw(24) << scientific << height.get(i) << " ";
	write_file << setw(24) << scientific << b.get(i) << " ";
	write_file << setw(24) << scientific << Qflow/pow(height.get(i),1.5) 
		   << " ";
	if (i!=(numpts+1))
	{ write_file << endl; }
    }
    write_file.close();

    // more output
    ofstream write_file2("bc.dat", ios::out);   
    write_file2 << setprecision(16);
    write_file2 << setw(24) << scientific << height.get(1) << endl;
    write_file2 << setw(24) << scientific << Qflow << endl;
    write_file2 << setw(24) << scientific << height.get(numpts+1) << endl;
    write_file2 << setw(24) << scientific << Qflow << endl;
    write_file2.close();

}

double bottom(double x)
{
    double b;

    b = 0.5e0*exp(-100.0e0*pow((x-0.5),2));
  
    return b;
}

double Newton(double a0,double a1,double a2,double guess,double TOL,int MaxIter)
{
    int i=0;
    int NumIter;
    int mstop=0;
    double root = guess;
    double fval,dfval;
    double err1,err2;

    fval = a0 + a1*guess + a2*pow(guess,2) + pow(guess,3);
  
    while (i<=MaxIter && mstop==0)    
    {
        dfval = a1 + 2.0e0*a2*guess + 3.0e0*pow(guess,2);
	root = guess - fval/dfval;

	err1 = fabs(root-guess);
	guess = root;
	fval = a0 + a1*guess + a2*pow(guess,2) + pow(guess,3);
	err2 = fabs(fval);

	if (err1<TOL && err2<TOL)
	{ mstop = 1; }

	i = i + 1;
    }
    
    if (mstop==0)
    {
        cout << endl;
	cout << "  ERROR: Newton Iteration did not converge ... " << endl;
	cout << "                      TOL = " << TOL << endl;
	cout << "                     err1 = " << err1 << endl;
	cout << "                     err2 = " << err2 << endl;
	cout << "     number of iterations = " << i << endl;
	cout << endl;
    }

    return root;
}
